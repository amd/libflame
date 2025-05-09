/* ./strsyl.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static logical c_false = FALSE_;
static integer c__2 = 2;
static real c_b26 = 1.f;
static real c_b30 = 0.f;
static logical c_true = TRUE_;
/* > \brief \b STRSYL */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download STRSYL + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strsyl.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strsyl.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strsyl.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE STRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, */
/* LDC, SCALE, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANA, TRANB */
/* INTEGER INFO, ISGN, LDA, LDB, LDC, M, N */
/* REAL SCALE */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), B( LDB, * ), C( LDC, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > STRSYL solves the real Sylvester matrix equation: */
/* > */
/* > op(A)*X + X*op(B) = scale*C or */
/* > op(A)*X - X*op(B) = scale*C, */
/* > */
/* > where op(A) = A or A**T, and A and B are both upper quasi- */
/* > triangular. A is M-by-M and B is N-by-N;
the right hand side C and */
/* > the solution X are M-by-N;
and scale is an output scale factor, set */
/* > <= 1 to avoid overflow in X. */
/* > */
/* > A and B must be in Schur canonical form (as returned by SHSEQR), that */
/* > is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
 */
/* > each 2-by-2 diagonal block has its diagonal elements equal and its */
/* > off-diagonal elements of opposite sign. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANA */
/* > \verbatim */
/* > TRANA is CHARACTER*1 */
/* > Specifies the option op(A): */
/* > = 'N': op(A) = A (No transpose) */
/* > = 'T': op(A) = A**T (Transpose) */
/* > = 'C': op(A) = A**H (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] TRANB */
/* > \verbatim */
/* > TRANB is CHARACTER*1 */
/* > Specifies the option op(B): */
/* > = 'N': op(B) = B (No transpose) */
/* > = 'T': op(B) = B**T (Transpose) */
/* > = 'C': op(B) = B**H (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] ISGN */
/* > \verbatim */
/* > ISGN is INTEGER */
/* > Specifies the sign in the equation: */
/* > = +1: solve op(A)*X + X*op(B) = scale*C */
/* > = -1: solve op(A)*X - X*op(B) = scale*C */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The order of the matrix A, and the number of rows in the */
/* > matrices X and C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix B, and the number of columns in the */
/* > matrices X and C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,M) */
/* > The upper quasi-triangular matrix A, in Schur canonical form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB,N) */
/* > The upper quasi-triangular matrix B, in Schur canonical form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is REAL array, dimension (LDC,N) */
/* > On entry, the M-by-N right hand side matrix C. */
/* > On exit, C is overwritten by the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= fla_max(1,M) */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* > SCALE is REAL */
/* > The scale factor, scale, set <= 1 to avoid overflow in X. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > = 1: A and B have common or very close eigenvalues;
perturbed */
/* > values were used to solve the equation (but the matrices */
/* > A and B are unchanged). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup trsyl */
/* ===================================================================== */
/* Subroutine */
void strsyl_(char *trana, char *tranb, integer *isgn, integer *m, integer *n, real *a, integer *lda,
             real *b, integer *ldb, real *c__, integer *ldc, real *scale, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF(
             "strsyl inputs: trana %c, tranb %c, isgn %" FLA_IS ", m %" FLA_IS ", n %" FLA_IS
             ", lda %" FLA_IS ", ldb %" FLA_IS ", ldc %" FLA_IS "",
             *trana, *tranb, *isgn, *m, *n, *lda, *ldb, *ldc);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;
    /* Local variables */
    integer j, k, l;
    real x[4] /* was [2][2] */
        ;
    integer k1, k2, l1, l2;
    real a11, db, da11, vec[4] /* was [2][2] */
        ,
        dum[1], eps, sgn;
    integer ierr;
    real smin;
    extern real sdot_(integer *, real *, integer *, real *, integer *);
    real suml, sumr;
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *);
    integer knext, lnext;
    real xnorm;
    extern /* Subroutine */
        void
        slaln2_(logical *, integer *, integer *, real *, real *, real *, integer *, real *, real *,
                real *, integer *, real *, real *, real *, integer *, real *, real *, integer *),
        slasy2_(logical *, logical *, integer *, integer *, integer *, real *, integer *, real *,
                integer *, real *, integer *, real *, real *, integer *, real *, integer *);
    real scaloc;
    extern real slamch_(char *), slange_(char *, integer *, integer *, real *, integer *, real *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real bignum;
    logical notrna, notrnb;
    real smlnum;
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
    /* Decode and Test input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    /* Function Body */
    notrna = lsame_(trana, "N", 1, 1);
    notrnb = lsame_(tranb, "N", 1, 1);
    *info = 0;
    if(!notrna && !lsame_(trana, "T", 1, 1) && !lsame_(trana, "C", 1, 1))
    {
        *info = -1;
    }
    else if(!notrnb && !lsame_(tranb, "T", 1, 1) && !lsame_(tranb, "C", 1, 1))
    {
        *info = -2;
    }
    else if(*isgn != 1 && *isgn != -1)
    {
        *info = -3;
    }
    else if(*m < 0)
    {
        *info = -4;
    }
    else if(*n < 0)
    {
        *info = -5;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -7;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -9;
    }
    else if(*ldc < fla_max(1, *m))
    {
        *info = -11;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("STRSYL", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    *scale = 1.f;
    if(*m == 0 || *n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Set constants to control overflow */
    eps = slamch_("P");
    smlnum = slamch_("S");
    bignum = 1.f / smlnum;
    smlnum = smlnum * (real)(*m * *n) / eps;
    bignum = 1.f / smlnum;
    /* Computing MAX */
    r__1 = smlnum, r__2 = eps * slange_("M", m, m, &a[a_offset], lda, dum);
    r__1 = fla_max(r__1, r__2);
    r__2 = eps * slange_("M", n, n, &b[b_offset], ldb, dum); // ; expr subst
    smin = fla_max(r__1, r__2);
    sgn = (real)(*isgn);
    if(notrna && notrnb)
    {
        /* Solve A*X + ISGN*X*B = scale*C. */
        /* The (K,L)th block of X is determined starting from */
        /* bottom-left corner column by column by */
        /* A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */
        /* Where */
        /* M L-1 */
        /* R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]. */
        /* I=K+1 J=1 */
        /* Start column loop (index = L) */
        /* L1 (L2) : column index of the first (first) row of X(K,L). */
        lnext = 1;
        i__1 = *n;
        for(l = 1; l <= i__1; ++l)
        {
            if(l < lnext)
            {
                goto L70;
            }
            if(l == *n)
            {
                l1 = l;
                l2 = l;
            }
            else
            {
                if(b[l + 1 + l * b_dim1] != 0.f)
                {
                    l1 = l;
                    l2 = l + 1;
                    lnext = l + 2;
                }
                else
                {
                    l1 = l;
                    l2 = l;
                    lnext = l + 1;
                }
            }
            /* Start row loop (index = K) */
            /* K1 (K2): row index of the first (last) row of X(K,L). */
            knext = *m;
            for(k = *m; k >= 1; --k)
            {
                if(k > knext)
                {
                    goto L60;
                }
                if(k == 1)
                {
                    k1 = k;
                    k2 = k;
                }
                else
                {
                    if(a[k + (k - 1) * a_dim1] != 0.f)
                    {
                        k1 = k - 1;
                        k2 = k;
                        knext = k - 2;
                    }
                    else
                    {
                        k1 = k;
                        k2 = k;
                        knext = k - 1;
                    }
                }
                if(l1 == l2 && k1 == k2)
                {
                    i__2 = *m - k1;
                    /* Computing MIN */
                    i__3 = k1 + 1;
                    /* Computing MIN */
                    i__4 = k1 + 1;
                    suml = sdot_(&i__2, &a[k1 + fla_min(i__3, *m) * a_dim1], lda,
                                 &c__[fla_min(i__4, *m) + l1 * c_dim1], &c__1);
                    i__2 = l1 - 1;
                    sumr = sdot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
                    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
                    scaloc = 1.f;
                    a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
                    da11 = f2c_abs(a11);
                    if(da11 <= smin)
                    {
                        a11 = smin;
                        da11 = smin;
                        *info = 1;
                    }
                    db = f2c_abs(vec[0]);
                    if(da11 < 1.f && db > 1.f)
                    {
                        if(db > bignum * da11)
                        {
                            scaloc = 1.f / db;
                        }
                    }
                    x[0] = vec[0] * scaloc / a11;
                    if(scaloc != 1.f)
                    {
                        i__2 = *n;
                        for(j = 1; j <= i__2; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L10: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                }
                else if(l1 == l2 && k1 != k2)
                {
                    i__2 = *m - k2;
                    /* Computing MIN */
                    i__3 = k2 + 1;
                    /* Computing MIN */
                    i__4 = k2 + 1;
                    suml = sdot_(&i__2, &a[k1 + fla_min(i__3, *m) * a_dim1], lda,
                                 &c__[fla_min(i__4, *m) + l1 * c_dim1], &c__1);
                    i__2 = l1 - 1;
                    sumr = sdot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
                    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
                    i__2 = *m - k2;
                    /* Computing MIN */
                    i__3 = k2 + 1;
                    /* Computing MIN */
                    i__4 = k2 + 1;
                    suml = sdot_(&i__2, &a[k2 + fla_min(i__3, *m) * a_dim1], lda,
                                 &c__[fla_min(i__4, *m) + l1 * c_dim1], &c__1);
                    i__2 = l1 - 1;
                    sumr = sdot_(&i__2, &c__[k2 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
                    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);
                    r__1 = -sgn * b[l1 + l1 * b_dim1];
                    slaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &a[k1 + k1 * a_dim1], lda,
                            &c_b26, &c_b26, vec, &c__2, &r__1, &c_b30, x, &c__2, &scaloc, &xnorm,
                            &ierr);
                    if(ierr != 0)
                    {
                        *info = 1;
                    }
                    if(scaloc != 1.f)
                    {
                        i__2 = *n;
                        for(j = 1; j <= i__2; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L20: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                    c__[k2 + l1 * c_dim1] = x[1];
                }
                else if(l1 != l2 && k1 == k2)
                {
                    i__2 = *m - k1;
                    /* Computing MIN */
                    i__3 = k1 + 1;
                    /* Computing MIN */
                    i__4 = k1 + 1;
                    suml = sdot_(&i__2, &a[k1 + fla_min(i__3, *m) * a_dim1], lda,
                                 &c__[fla_min(i__4, *m) + l1 * c_dim1], &c__1);
                    i__2 = l1 - 1;
                    sumr = sdot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
                    vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * sumr));
                    i__2 = *m - k1;
                    /* Computing MIN */
                    i__3 = k1 + 1;
                    /* Computing MIN */
                    i__4 = k1 + 1;
                    suml = sdot_(&i__2, &a[k1 + fla_min(i__3, *m) * a_dim1], lda,
                                 &c__[fla_min(i__4, *m) + l2 * c_dim1], &c__1);
                    i__2 = l1 - 1;
                    sumr = sdot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l2 * b_dim1 + 1], &c__1);
                    vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * sumr));
                    r__1 = -sgn * a[k1 + k1 * a_dim1];
                    slaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &b[l1 + l1 * b_dim1], ldb, &c_b26,
                            &c_b26, vec, &c__2, &r__1, &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
                    if(ierr != 0)
                    {
                        *info = 1;
                    }
                    if(scaloc != 1.f)
                    {
                        i__2 = *n;
                        for(j = 1; j <= i__2; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L40: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                    c__[k1 + l2 * c_dim1] = x[1];
                }
                else if(l1 != l2 && k1 != k2)
                {
                    i__2 = *m - k2;
                    /* Computing MIN */
                    i__3 = k2 + 1;
                    /* Computing MIN */
                    i__4 = k2 + 1;
                    suml = sdot_(&i__2, &a[k1 + fla_min(i__3, *m) * a_dim1], lda,
                                 &c__[fla_min(i__4, *m) + l1 * c_dim1], &c__1);
                    i__2 = l1 - 1;
                    sumr = sdot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
                    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
                    i__2 = *m - k2;
                    /* Computing MIN */
                    i__3 = k2 + 1;
                    /* Computing MIN */
                    i__4 = k2 + 1;
                    suml = sdot_(&i__2, &a[k1 + fla_min(i__3, *m) * a_dim1], lda,
                                 &c__[fla_min(i__4, *m) + l2 * c_dim1], &c__1);
                    i__2 = l1 - 1;
                    sumr = sdot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l2 * b_dim1 + 1], &c__1);
                    vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);
                    i__2 = *m - k2;
                    /* Computing MIN */
                    i__3 = k2 + 1;
                    /* Computing MIN */
                    i__4 = k2 + 1;
                    suml = sdot_(&i__2, &a[k2 + fla_min(i__3, *m) * a_dim1], lda,
                                 &c__[fla_min(i__4, *m) + l1 * c_dim1], &c__1);
                    i__2 = l1 - 1;
                    sumr = sdot_(&i__2, &c__[k2 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
                    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);
                    i__2 = *m - k2;
                    /* Computing MIN */
                    i__3 = k2 + 1;
                    /* Computing MIN */
                    i__4 = k2 + 1;
                    suml = sdot_(&i__2, &a[k2 + fla_min(i__3, *m) * a_dim1], lda,
                                 &c__[fla_min(i__4, *m) + l2 * c_dim1], &c__1);
                    i__2 = l1 - 1;
                    sumr = sdot_(&i__2, &c__[k2 + c_dim1], ldc, &b[l2 * b_dim1 + 1], &c__1);
                    vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);
                    slasy2_(&c_false, &c_false, isgn, &c__2, &c__2, &a[k1 + k1 * a_dim1], lda,
                            &b[l1 + l1 * b_dim1], ldb, vec, &c__2, &scaloc, x, &c__2, &xnorm,
                            &ierr);
                    if(ierr != 0)
                    {
                        *info = 1;
                    }
                    if(scaloc != 1.f)
                    {
                        i__2 = *n;
                        for(j = 1; j <= i__2; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L50: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                    c__[k1 + l2 * c_dim1] = x[2];
                    c__[k2 + l1 * c_dim1] = x[1];
                    c__[k2 + l2 * c_dim1] = x[3];
                }
            L60:;
            }
        L70:;
        }
    }
    else if(!notrna && notrnb)
    {
        /* Solve A**T *X + ISGN*X*B = scale*C. */
        /* The (K,L)th block of X is determined starting from */
        /* upper-left corner column by column by */
        /* A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */
        /* Where */
        /* K-1 L-1 */
        /* R(K,L) = SUM [A(I,K)**T*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)] */
        /* I=1 J=1 */
        /* Start column loop (index = L) */
        /* L1 (L2): column index of the first (last) row of X(K,L) */
        lnext = 1;
        i__1 = *n;
        for(l = 1; l <= i__1; ++l)
        {
            if(l < lnext)
            {
                goto L130;
            }
            if(l == *n)
            {
                l1 = l;
                l2 = l;
            }
            else
            {
                if(b[l + 1 + l * b_dim1] != 0.f)
                {
                    l1 = l;
                    l2 = l + 1;
                    lnext = l + 2;
                }
                else
                {
                    l1 = l;
                    l2 = l;
                    lnext = l + 1;
                }
            }
            /* Start row loop (index = K) */
            /* K1 (K2): row index of the first (last) row of X(K,L) */
            knext = 1;
            i__2 = *m;
            for(k = 1; k <= i__2; ++k)
            {
                if(k < knext)
                {
                    goto L120;
                }
                if(k == *m)
                {
                    k1 = k;
                    k2 = k;
                }
                else
                {
                    if(a[k + 1 + k * a_dim1] != 0.f)
                    {
                        k1 = k;
                        k2 = k + 1;
                        knext = k + 2;
                    }
                    else
                    {
                        k1 = k;
                        k2 = k;
                        knext = k + 1;
                    }
                }
                if(l1 == l2 && k1 == k2)
                {
                    i__3 = k1 - 1;
                    suml = sdot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
                    i__3 = l1 - 1;
                    sumr = sdot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
                    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
                    scaloc = 1.f;
                    a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
                    da11 = f2c_abs(a11);
                    if(da11 <= smin)
                    {
                        a11 = smin;
                        da11 = smin;
                        *info = 1;
                    }
                    db = f2c_abs(vec[0]);
                    if(da11 < 1.f && db > 1.f)
                    {
                        if(db > bignum * da11)
                        {
                            scaloc = 1.f / db;
                        }
                    }
                    x[0] = vec[0] * scaloc / a11;
                    if(scaloc != 1.f)
                    {
                        i__3 = *n;
                        for(j = 1; j <= i__3; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L80: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                }
                else if(l1 == l2 && k1 != k2)
                {
                    i__3 = k1 - 1;
                    suml = sdot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
                    i__3 = l1 - 1;
                    sumr = sdot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
                    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
                    i__3 = k1 - 1;
                    suml = sdot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
                    i__3 = l1 - 1;
                    sumr = sdot_(&i__3, &c__[k2 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
                    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);
                    r__1 = -sgn * b[l1 + l1 * b_dim1];
                    slaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &a[k1 + k1 * a_dim1], lda, &c_b26,
                            &c_b26, vec, &c__2, &r__1, &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
                    if(ierr != 0)
                    {
                        *info = 1;
                    }
                    if(scaloc != 1.f)
                    {
                        i__3 = *n;
                        for(j = 1; j <= i__3; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L90: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                    c__[k2 + l1 * c_dim1] = x[1];
                }
                else if(l1 != l2 && k1 == k2)
                {
                    i__3 = k1 - 1;
                    suml = sdot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
                    i__3 = l1 - 1;
                    sumr = sdot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
                    vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * sumr));
                    i__3 = k1 - 1;
                    suml = sdot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1);
                    i__3 = l1 - 1;
                    sumr = sdot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l2 * b_dim1 + 1], &c__1);
                    vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * sumr));
                    r__1 = -sgn * a[k1 + k1 * a_dim1];
                    slaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &b[l1 + l1 * b_dim1], ldb, &c_b26,
                            &c_b26, vec, &c__2, &r__1, &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
                    if(ierr != 0)
                    {
                        *info = 1;
                    }
                    if(scaloc != 1.f)
                    {
                        i__3 = *n;
                        for(j = 1; j <= i__3; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L100: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                    c__[k1 + l2 * c_dim1] = x[1];
                }
                else if(l1 != l2 && k1 != k2)
                {
                    i__3 = k1 - 1;
                    suml = sdot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
                    i__3 = l1 - 1;
                    sumr = sdot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
                    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
                    i__3 = k1 - 1;
                    suml = sdot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1);
                    i__3 = l1 - 1;
                    sumr = sdot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l2 * b_dim1 + 1], &c__1);
                    vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);
                    i__3 = k1 - 1;
                    suml = sdot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
                    i__3 = l1 - 1;
                    sumr = sdot_(&i__3, &c__[k2 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
                    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);
                    i__3 = k1 - 1;
                    suml = sdot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1);
                    i__3 = l1 - 1;
                    sumr = sdot_(&i__3, &c__[k2 + c_dim1], ldc, &b[l2 * b_dim1 + 1], &c__1);
                    vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);
                    slasy2_(&c_true, &c_false, isgn, &c__2, &c__2, &a[k1 + k1 * a_dim1], lda,
                            &b[l1 + l1 * b_dim1], ldb, vec, &c__2, &scaloc, x, &c__2, &xnorm,
                            &ierr);
                    if(ierr != 0)
                    {
                        *info = 1;
                    }
                    if(scaloc != 1.f)
                    {
                        i__3 = *n;
                        for(j = 1; j <= i__3; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L110: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                    c__[k1 + l2 * c_dim1] = x[2];
                    c__[k2 + l1 * c_dim1] = x[1];
                    c__[k2 + l2 * c_dim1] = x[3];
                }
            L120:;
            }
        L130:;
        }
    }
    else if(!notrna && !notrnb)
    {
        /* Solve A**T*X + ISGN*X*B**T = scale*C. */
        /* The (K,L)th block of X is determined starting from */
        /* top-right corner column by column by */
        /* A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L) */
        /* Where */
        /* K-1 N */
        /* R(K,L) = SUM [A(I,K)**T*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T]. */
        /* I=1 J=L+1 */
        /* Start column loop (index = L) */
        /* L1 (L2): column index of the first (last) row of X(K,L) */
        lnext = *n;
        for(l = *n; l >= 1; --l)
        {
            if(l > lnext)
            {
                goto L190;
            }
            if(l == 1)
            {
                l1 = l;
                l2 = l;
            }
            else
            {
                if(b[l + (l - 1) * b_dim1] != 0.f)
                {
                    l1 = l - 1;
                    l2 = l;
                    lnext = l - 2;
                }
                else
                {
                    l1 = l;
                    l2 = l;
                    lnext = l - 1;
                }
            }
            /* Start row loop (index = K) */
            /* K1 (K2): row index of the first (last) row of X(K,L) */
            knext = 1;
            i__1 = *m;
            for(k = 1; k <= i__1; ++k)
            {
                if(k < knext)
                {
                    goto L180;
                }
                if(k == *m)
                {
                    k1 = k;
                    k2 = k;
                }
                else
                {
                    if(a[k + 1 + k * a_dim1] != 0.f)
                    {
                        k1 = k;
                        k2 = k + 1;
                        knext = k + 2;
                    }
                    else
                    {
                        k1 = k;
                        k2 = k;
                        knext = k + 1;
                    }
                }
                if(l1 == l2 && k1 == k2)
                {
                    i__2 = k1 - 1;
                    suml = sdot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
                    i__2 = *n - l1;
                    /* Computing MIN */
                    i__3 = l1 + 1;
                    /* Computing MIN */
                    i__4 = l1 + 1;
                    sumr = sdot_(&i__2, &c__[k1 + fla_min(i__3, *n) * c_dim1], ldc,
                                 &b[l1 + fla_min(i__4, *n) * b_dim1], ldb);
                    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
                    scaloc = 1.f;
                    a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
                    da11 = f2c_abs(a11);
                    if(da11 <= smin)
                    {
                        a11 = smin;
                        da11 = smin;
                        *info = 1;
                    }
                    db = f2c_abs(vec[0]);
                    if(da11 < 1.f && db > 1.f)
                    {
                        if(db > bignum * da11)
                        {
                            scaloc = 1.f / db;
                        }
                    }
                    x[0] = vec[0] * scaloc / a11;
                    if(scaloc != 1.f)
                    {
                        i__2 = *n;
                        for(j = 1; j <= i__2; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L140: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                }
                else if(l1 == l2 && k1 != k2)
                {
                    i__2 = k1 - 1;
                    suml = sdot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
                    i__2 = *n - l2;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    /* Computing MIN */
                    i__4 = l2 + 1;
                    sumr = sdot_(&i__2, &c__[k1 + fla_min(i__3, *n) * c_dim1], ldc,
                                 &b[l1 + fla_min(i__4, *n) * b_dim1], ldb);
                    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
                    i__2 = k1 - 1;
                    suml = sdot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
                    i__2 = *n - l2;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    /* Computing MIN */
                    i__4 = l2 + 1;
                    sumr = sdot_(&i__2, &c__[k2 + fla_min(i__3, *n) * c_dim1], ldc,
                                 &b[l1 + fla_min(i__4, *n) * b_dim1], ldb);
                    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);
                    r__1 = -sgn * b[l1 + l1 * b_dim1];
                    slaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &a[k1 + k1 * a_dim1], lda, &c_b26,
                            &c_b26, vec, &c__2, &r__1, &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
                    if(ierr != 0)
                    {
                        *info = 1;
                    }
                    if(scaloc != 1.f)
                    {
                        i__2 = *n;
                        for(j = 1; j <= i__2; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L150: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                    c__[k2 + l1 * c_dim1] = x[1];
                }
                else if(l1 != l2 && k1 == k2)
                {
                    i__2 = k1 - 1;
                    suml = sdot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
                    i__2 = *n - l2;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    /* Computing MIN */
                    i__4 = l2 + 1;
                    sumr = sdot_(&i__2, &c__[k1 + fla_min(i__3, *n) * c_dim1], ldc,
                                 &b[l1 + fla_min(i__4, *n) * b_dim1], ldb);
                    vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * sumr));
                    i__2 = k1 - 1;
                    suml = sdot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1);
                    i__2 = *n - l2;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    /* Computing MIN */
                    i__4 = l2 + 1;
                    sumr = sdot_(&i__2, &c__[k1 + fla_min(i__3, *n) * c_dim1], ldc,
                                 &b[l2 + fla_min(i__4, *n) * b_dim1], ldb);
                    vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * sumr));
                    r__1 = -sgn * a[k1 + k1 * a_dim1];
                    slaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &b[l1 + l1 * b_dim1], ldb,
                            &c_b26, &c_b26, vec, &c__2, &r__1, &c_b30, x, &c__2, &scaloc, &xnorm,
                            &ierr);
                    if(ierr != 0)
                    {
                        *info = 1;
                    }
                    if(scaloc != 1.f)
                    {
                        i__2 = *n;
                        for(j = 1; j <= i__2; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L160: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                    c__[k1 + l2 * c_dim1] = x[1];
                }
                else if(l1 != l2 && k1 != k2)
                {
                    i__2 = k1 - 1;
                    suml = sdot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
                    i__2 = *n - l2;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    /* Computing MIN */
                    i__4 = l2 + 1;
                    sumr = sdot_(&i__2, &c__[k1 + fla_min(i__3, *n) * c_dim1], ldc,
                                 &b[l1 + fla_min(i__4, *n) * b_dim1], ldb);
                    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
                    i__2 = k1 - 1;
                    suml = sdot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1);
                    i__2 = *n - l2;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    /* Computing MIN */
                    i__4 = l2 + 1;
                    sumr = sdot_(&i__2, &c__[k1 + fla_min(i__3, *n) * c_dim1], ldc,
                                 &b[l2 + fla_min(i__4, *n) * b_dim1], ldb);
                    vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);
                    i__2 = k1 - 1;
                    suml = sdot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
                    i__2 = *n - l2;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    /* Computing MIN */
                    i__4 = l2 + 1;
                    sumr = sdot_(&i__2, &c__[k2 + fla_min(i__3, *n) * c_dim1], ldc,
                                 &b[l1 + fla_min(i__4, *n) * b_dim1], ldb);
                    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);
                    i__2 = k1 - 1;
                    suml = sdot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1);
                    i__2 = *n - l2;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    /* Computing MIN */
                    i__4 = l2 + 1;
                    sumr = sdot_(&i__2, &c__[k2 + fla_min(i__3, *n) * c_dim1], ldc,
                                 &b[l2 + fla_min(i__4, *n) * b_dim1], ldb);
                    vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);
                    slasy2_(&c_true, &c_true, isgn, &c__2, &c__2, &a[k1 + k1 * a_dim1], lda,
                            &b[l1 + l1 * b_dim1], ldb, vec, &c__2, &scaloc, x, &c__2, &xnorm,
                            &ierr);
                    if(ierr != 0)
                    {
                        *info = 1;
                    }
                    if(scaloc != 1.f)
                    {
                        i__2 = *n;
                        for(j = 1; j <= i__2; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L170: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                    c__[k1 + l2 * c_dim1] = x[2];
                    c__[k2 + l1 * c_dim1] = x[1];
                    c__[k2 + l2 * c_dim1] = x[3];
                }
            L180:;
            }
        L190:;
        }
    }
    else if(notrna && !notrnb)
    {
        /* Solve A*X + ISGN*X*B**T = scale*C. */
        /* The (K,L)th block of X is determined starting from */
        /* bottom-right corner column by column by */
        /* A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L) */
        /* Where */
        /* M N */
        /* R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T]. */
        /* I=K+1 J=L+1 */
        /* Start column loop (index = L) */
        /* L1 (L2): column index of the first (last) row of X(K,L) */
        lnext = *n;
        for(l = *n; l >= 1; --l)
        {
            if(l > lnext)
            {
                goto L250;
            }
            if(l == 1)
            {
                l1 = l;
                l2 = l;
            }
            else
            {
                if(b[l + (l - 1) * b_dim1] != 0.f)
                {
                    l1 = l - 1;
                    l2 = l;
                    lnext = l - 2;
                }
                else
                {
                    l1 = l;
                    l2 = l;
                    lnext = l - 1;
                }
            }
            /* Start row loop (index = K) */
            /* K1 (K2): row index of the first (last) row of X(K,L) */
            knext = *m;
            for(k = *m; k >= 1; --k)
            {
                if(k > knext)
                {
                    goto L240;
                }
                if(k == 1)
                {
                    k1 = k;
                    k2 = k;
                }
                else
                {
                    if(a[k + (k - 1) * a_dim1] != 0.f)
                    {
                        k1 = k - 1;
                        k2 = k;
                        knext = k - 2;
                    }
                    else
                    {
                        k1 = k;
                        k2 = k;
                        knext = k - 1;
                    }
                }
                if(l1 == l2 && k1 == k2)
                {
                    i__1 = *m - k1;
                    /* Computing MIN */
                    i__2 = k1 + 1;
                    /* Computing MIN */
                    i__3 = k1 + 1;
                    suml = sdot_(&i__1, &a[k1 + fla_min(i__2, *m) * a_dim1], lda,
                                 &c__[fla_min(i__3, *m) + l1 * c_dim1], &c__1);
                    i__1 = *n - l1;
                    /* Computing MIN */
                    i__2 = l1 + 1;
                    /* Computing MIN */
                    i__3 = l1 + 1;
                    sumr = sdot_(&i__1, &c__[k1 + fla_min(i__2, *n) * c_dim1], ldc,
                                 &b[l1 + fla_min(i__3, *n) * b_dim1], ldb);
                    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
                    scaloc = 1.f;
                    a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
                    da11 = f2c_abs(a11);
                    if(da11 <= smin)
                    {
                        a11 = smin;
                        da11 = smin;
                        *info = 1;
                    }
                    db = f2c_abs(vec[0]);
                    if(da11 < 1.f && db > 1.f)
                    {
                        if(db > bignum * da11)
                        {
                            scaloc = 1.f / db;
                        }
                    }
                    x[0] = vec[0] * scaloc / a11;
                    if(scaloc != 1.f)
                    {
                        i__1 = *n;
                        for(j = 1; j <= i__1; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L200: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                }
                else if(l1 == l2 && k1 != k2)
                {
                    i__1 = *m - k2;
                    /* Computing MIN */
                    i__2 = k2 + 1;
                    /* Computing MIN */
                    i__3 = k2 + 1;
                    suml = sdot_(&i__1, &a[k1 + fla_min(i__2, *m) * a_dim1], lda,
                                 &c__[fla_min(i__3, *m) + l1 * c_dim1], &c__1);
                    i__1 = *n - l2;
                    /* Computing MIN */
                    i__2 = l2 + 1;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    sumr = sdot_(&i__1, &c__[k1 + fla_min(i__2, *n) * c_dim1], ldc,
                                 &b[l1 + fla_min(i__3, *n) * b_dim1], ldb);
                    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
                    i__1 = *m - k2;
                    /* Computing MIN */
                    i__2 = k2 + 1;
                    /* Computing MIN */
                    i__3 = k2 + 1;
                    suml = sdot_(&i__1, &a[k2 + fla_min(i__2, *m) * a_dim1], lda,
                                 &c__[fla_min(i__3, *m) + l1 * c_dim1], &c__1);
                    i__1 = *n - l2;
                    /* Computing MIN */
                    i__2 = l2 + 1;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    sumr = sdot_(&i__1, &c__[k2 + fla_min(i__2, *n) * c_dim1], ldc,
                                 &b[l1 + fla_min(i__3, *n) * b_dim1], ldb);
                    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);
                    r__1 = -sgn * b[l1 + l1 * b_dim1];
                    slaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &a[k1 + k1 * a_dim1], lda,
                            &c_b26, &c_b26, vec, &c__2, &r__1, &c_b30, x, &c__2, &scaloc, &xnorm,
                            &ierr);
                    if(ierr != 0)
                    {
                        *info = 1;
                    }
                    if(scaloc != 1.f)
                    {
                        i__1 = *n;
                        for(j = 1; j <= i__1; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L210: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                    c__[k2 + l1 * c_dim1] = x[1];
                }
                else if(l1 != l2 && k1 == k2)
                {
                    i__1 = *m - k1;
                    /* Computing MIN */
                    i__2 = k1 + 1;
                    /* Computing MIN */
                    i__3 = k1 + 1;
                    suml = sdot_(&i__1, &a[k1 + fla_min(i__2, *m) * a_dim1], lda,
                                 &c__[fla_min(i__3, *m) + l1 * c_dim1], &c__1);
                    i__1 = *n - l2;
                    /* Computing MIN */
                    i__2 = l2 + 1;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    sumr = sdot_(&i__1, &c__[k1 + fla_min(i__2, *n) * c_dim1], ldc,
                                 &b[l1 + fla_min(i__3, *n) * b_dim1], ldb);
                    vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * sumr));
                    i__1 = *m - k1;
                    /* Computing MIN */
                    i__2 = k1 + 1;
                    /* Computing MIN */
                    i__3 = k1 + 1;
                    suml = sdot_(&i__1, &a[k1 + fla_min(i__2, *m) * a_dim1], lda,
                                 &c__[fla_min(i__3, *m) + l2 * c_dim1], &c__1);
                    i__1 = *n - l2;
                    /* Computing MIN */
                    i__2 = l2 + 1;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    sumr = sdot_(&i__1, &c__[k1 + fla_min(i__2, *n) * c_dim1], ldc,
                                 &b[l2 + fla_min(i__3, *n) * b_dim1], ldb);
                    vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * sumr));
                    r__1 = -sgn * a[k1 + k1 * a_dim1];
                    slaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &b[l1 + l1 * b_dim1], ldb,
                            &c_b26, &c_b26, vec, &c__2, &r__1, &c_b30, x, &c__2, &scaloc, &xnorm,
                            &ierr);
                    if(ierr != 0)
                    {
                        *info = 1;
                    }
                    if(scaloc != 1.f)
                    {
                        i__1 = *n;
                        for(j = 1; j <= i__1; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L220: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                    c__[k1 + l2 * c_dim1] = x[1];
                }
                else if(l1 != l2 && k1 != k2)
                {
                    i__1 = *m - k2;
                    /* Computing MIN */
                    i__2 = k2 + 1;
                    /* Computing MIN */
                    i__3 = k2 + 1;
                    suml = sdot_(&i__1, &a[k1 + fla_min(i__2, *m) * a_dim1], lda,
                                 &c__[fla_min(i__3, *m) + l1 * c_dim1], &c__1);
                    i__1 = *n - l2;
                    /* Computing MIN */
                    i__2 = l2 + 1;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    sumr = sdot_(&i__1, &c__[k1 + fla_min(i__2, *n) * c_dim1], ldc,
                                 &b[l1 + fla_min(i__3, *n) * b_dim1], ldb);
                    vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
                    i__1 = *m - k2;
                    /* Computing MIN */
                    i__2 = k2 + 1;
                    /* Computing MIN */
                    i__3 = k2 + 1;
                    suml = sdot_(&i__1, &a[k1 + fla_min(i__2, *m) * a_dim1], lda,
                                 &c__[fla_min(i__3, *m) + l2 * c_dim1], &c__1);
                    i__1 = *n - l2;
                    /* Computing MIN */
                    i__2 = l2 + 1;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    sumr = sdot_(&i__1, &c__[k1 + fla_min(i__2, *n) * c_dim1], ldc,
                                 &b[l2 + fla_min(i__3, *n) * b_dim1], ldb);
                    vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);
                    i__1 = *m - k2;
                    /* Computing MIN */
                    i__2 = k2 + 1;
                    /* Computing MIN */
                    i__3 = k2 + 1;
                    suml = sdot_(&i__1, &a[k2 + fla_min(i__2, *m) * a_dim1], lda,
                                 &c__[fla_min(i__3, *m) + l1 * c_dim1], &c__1);
                    i__1 = *n - l2;
                    /* Computing MIN */
                    i__2 = l2 + 1;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    sumr = sdot_(&i__1, &c__[k2 + fla_min(i__2, *n) * c_dim1], ldc,
                                 &b[l1 + fla_min(i__3, *n) * b_dim1], ldb);
                    vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);
                    i__1 = *m - k2;
                    /* Computing MIN */
                    i__2 = k2 + 1;
                    /* Computing MIN */
                    i__3 = k2 + 1;
                    suml = sdot_(&i__1, &a[k2 + fla_min(i__2, *m) * a_dim1], lda,
                                 &c__[fla_min(i__3, *m) + l2 * c_dim1], &c__1);
                    i__1 = *n - l2;
                    /* Computing MIN */
                    i__2 = l2 + 1;
                    /* Computing MIN */
                    i__3 = l2 + 1;
                    sumr = sdot_(&i__1, &c__[k2 + fla_min(i__2, *n) * c_dim1], ldc,
                                 &b[l2 + fla_min(i__3, *n) * b_dim1], ldb);
                    vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);
                    slasy2_(&c_false, &c_true, isgn, &c__2, &c__2, &a[k1 + k1 * a_dim1], lda,
                            &b[l1 + l1 * b_dim1], ldb, vec, &c__2, &scaloc, x, &c__2, &xnorm,
                            &ierr);
                    if(ierr != 0)
                    {
                        *info = 1;
                    }
                    if(scaloc != 1.f)
                    {
                        i__1 = *n;
                        for(j = 1; j <= i__1; ++j)
                        {
                            sscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                            /* L230: */
                        }
                        *scale *= scaloc;
                    }
                    c__[k1 + l1 * c_dim1] = x[0];
                    c__[k1 + l2 * c_dim1] = x[2];
                    c__[k2 + l1 * c_dim1] = x[1];
                    c__[k2 + l2 * c_dim1] = x[3];
                }
            L240:;
            }
        L250:;
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of STRSYL */
}
/* strsyl_ */
