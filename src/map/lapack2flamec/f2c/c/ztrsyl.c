/* ./ztrsyl.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief \b ZTRSYL */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZTRSYL + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrsyl.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrsyl.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrsyl.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, */
/* LDC, SCALE, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANA, TRANB */
/* INTEGER INFO, ISGN, LDA, LDB, LDC, M, N */
/* DOUBLE PRECISION SCALE */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), B( LDB, * ), C( LDC, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTRSYL solves the scomplex Sylvester matrix equation: */
/* > */
/* > op(A)*X + X*op(B) = scale*C or */
/* > op(A)*X - X*op(B) = scale*C, */
/* > */
/* > where op(A) = A or A**H, and A and B are both upper triangular. A is */
/* > M-by-M and B is N-by-N;
the right hand side C and the solution X are */
/* > M-by-N;
and scale is an output scale factor, set <= 1 to avoid */
/* > overflow in X. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANA */
/* > \verbatim */
/* > TRANA is CHARACTER*1 */
/* > Specifies the option op(A): */
/* > = 'N': op(A) = A (No transpose) */
/* > = 'C': op(A) = A**H (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] TRANB */
/* > \verbatim */
/* > TRANB is CHARACTER*1 */
/* > Specifies the option op(B): */
/* > = 'N': op(B) = B (No transpose) */
/* > = 'C': op(B) = B**H (Conjugate transpose) */
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
/* > A is COMPLEX*16 array, dimension (LDA,M) */
/* > The upper triangular matrix A. */
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
/* > B is COMPLEX*16 array, dimension (LDB,N) */
/* > The upper triangular matrix B. */
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
/* > C is COMPLEX*16 array, dimension (LDC,N) */
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
/* > SCALE is DOUBLE PRECISION */
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
/** Generated wrapper function */
void ztrsyl_(char *trana, char *tranb, aocl_int_t *isgn, aocl_int_t *m, aocl_int_t *n,
             dcomplex *a, aocl_int_t *lda, dcomplex *b, aocl_int_t *ldb,
             dcomplex *c__, aocl_int_t *ldc, doublereal *scale, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_ztrsyl(trana, tranb, isgn, m, n, a, lda, b, ldb, c__, ldc, scale, info);
#else
    aocl_int64_t isgn_64 = *isgn;
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t ldc_64 = *ldc;
    aocl_int64_t info_64 = *info;

    aocl_lapack_ztrsyl(trana, tranb, &isgn_64, &m_64, &n_64, a, &lda_64, b, &ldb_64, c__, &ldc_64,
                       scale, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_ztrsyl(char *trana, char *tranb, aocl_int64_t *isgn, aocl_int64_t *m,
                        aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                        aocl_int64_t *ldb, dcomplex *c__, aocl_int64_t *ldc, doublereal *scale,
                        aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("ztrsyl inputs: trana %c, tranb %c, isgn %" FLA_IS ", m %" FLA_IS
                      ", n %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS ", ldc %" FLA_IS "",
                      *trana, *tranb, *isgn, *m, *n, *lda, *ldb, *ldc);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    dcomplex z__1, z__2, z__3, z__4;
    /* Builtin functions */
    double d_imag(dcomplex *);
    void d_cnjg(dcomplex *, dcomplex *);
    /* Local variables */
    aocl_int64_t j, k, l;
    dcomplex a11;
    doublereal db;
    dcomplex x11;
    doublereal da11;
    dcomplex vec;
    doublereal dum[1], eps, sgn, smin;
    dcomplex suml, sumr;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    extern doublereal dlamch_(char *);
    doublereal scaloc;
    doublereal bignum;
    extern /* Double Complex */
        void
        zladiv_f2c_(dcomplex *, dcomplex *, dcomplex *);
    logical notrna, notrnb;
    doublereal smlnum;
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
    if(!notrna && !lsame_(trana, "C", 1, 1))
    {
        *info = -1;
    }
    else if(!notrnb && !lsame_(tranb, "C", 1, 1))
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
        aocl_blas_xerbla("ZTRSYL", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    *scale = 1.;
    if(*m == 0 || *n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Set constants to control overflow */
    eps = dlamch_("P");
    smlnum = dlamch_("S");
    bignum = 1. / smlnum;
    smlnum = smlnum * (doublereal)(*m * *n) / eps;
    bignum = 1. / smlnum;
    /* Computing MAX */
    d__1 = smlnum, d__2 = eps * aocl_lapack_zlange("M", m, m, &a[a_offset], lda, dum);
    d__1 = fla_max(d__1, d__2);
    d__2 = eps * aocl_lapack_zlange("M", n, n, &b[b_offset], ldb, dum); // ; expr subst
    smin = fla_max(d__1, d__2);
    sgn = (doublereal)(*isgn);
    if(notrna && notrnb)
    {
        /* Solve A*X + ISGN*X*B = scale*C. */
        /* The (K,L)th block of X is determined starting from */
        /* bottom-left corner column by column by */
        /* A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */
        /* Where */
        /* M L-1 */
        /* R(K,L) = SUM [A(K,I)*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]. */
        /* I=K+1 J=1 */
        i__1 = *n;
        for(l = 1; l <= i__1; ++l)
        {
            for(k = *m; k >= 1; --k)
            {
                i__2 = *m - k;
                /* Computing MIN */
                i__3 = k + 1;
                /* Computing MIN */
                i__4 = k + 1;
                aocl_lapack_zdotu_f2c(&z__1, &i__2, &a[k + fla_min(i__3, *m) * a_dim1], lda,
                           &c__[fla_min(i__4, *m) + l * c_dim1], &c__1);
                suml.real = z__1.real;
                suml.imag = z__1.imag; // , expr subst
                i__2 = l - 1;
                aocl_lapack_zdotu_f2c(&z__1, &i__2, &c__[k + c_dim1], ldc, &b[l * b_dim1 + 1], &c__1);
                sumr.real = z__1.real;
                sumr.imag = z__1.imag; // , expr subst
                i__2 = k + l * c_dim1;
                z__3.real = sgn * sumr.real;
                z__3.imag = sgn * sumr.imag; // , expr subst
                z__2.real = suml.real + z__3.real;
                z__2.imag = suml.imag + z__3.imag; // , expr subst
                z__1.real = c__[i__2].real - z__2.real;
                z__1.imag = c__[i__2].imag - z__2.imag; // , expr subst
                vec.real = z__1.real;
                vec.imag = z__1.imag; // , expr subst
                scaloc = 1.;
                i__2 = k + k * a_dim1;
                i__3 = l + l * b_dim1;
                z__2.real = sgn * b[i__3].real;
                z__2.imag = sgn * b[i__3].imag; // , expr subst
                z__1.real = a[i__2].real + z__2.real;
                z__1.imag = a[i__2].imag + z__2.imag; // , expr subst
                a11.real = z__1.real;
                a11.imag = z__1.imag; // , expr subst
                da11 = (d__1 = a11.real, f2c_dabs(d__1)) + (d__2 = d_imag(&a11), f2c_dabs(d__2));
                if(da11 <= smin)
                {
                    a11.real = smin;
                    a11.imag = 0.; // , expr subst
                    da11 = smin;
                    *info = 1;
                }
                db = (d__1 = vec.real, f2c_dabs(d__1)) + (d__2 = d_imag(&vec), f2c_dabs(d__2));
                if(da11 < 1. && db > 1.)
                {
                    if(db > bignum * da11)
                    {
                        scaloc = 1. / db;
                    }
                }
                z__3.real = scaloc;
                z__3.imag = 0.; // , expr subst
                z__2.real = vec.real * z__3.real - vec.imag * z__3.imag;
                z__2.imag = vec.real * z__3.imag + vec.imag * z__3.real; // , expr subst
                zladiv_f2c_(&z__1, &z__2, &a11);
                x11.real = z__1.real;
                x11.imag = z__1.imag; // , expr subst
                if(scaloc != 1.)
                {
                    i__2 = *n;
                    for(j = 1; j <= i__2; ++j)
                    {
                        aocl_blas_zdscal(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                        /* L10: */
                    }
                    *scale *= scaloc;
                }
                i__2 = k + l * c_dim1;
                c__[i__2].real = x11.real;
                c__[i__2].imag = x11.imag; // , expr subst
                /* L20: */
            }
            /* L30: */
        }
    }
    else if(!notrna && notrnb)
    {
        /* Solve A**H *X + ISGN*X*B = scale*C. */
        /* The (K,L)th block of X is determined starting from */
        /* upper-left corner column by column by */
        /* A**H(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */
        /* Where */
        /* K-1 L-1 */
        /* R(K,L) = SUM [A**H(I,K)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)] */
        /* I=1 J=1 */
        i__1 = *n;
        for(l = 1; l <= i__1; ++l)
        {
            i__2 = *m;
            for(k = 1; k <= i__2; ++k)
            {
                i__3 = k - 1;
                aocl_lapack_zdotc_f2c(&z__1, &i__3, &a[k * a_dim1 + 1], &c__1, &c__[l * c_dim1 + 1], &c__1);
                suml.real = z__1.real;
                suml.imag = z__1.imag; // , expr subst
                i__3 = l - 1;
                aocl_lapack_zdotu_f2c(&z__1, &i__3, &c__[k + c_dim1], ldc, &b[l * b_dim1 + 1], &c__1);
                sumr.real = z__1.real;
                sumr.imag = z__1.imag; // , expr subst
                i__3 = k + l * c_dim1;
                z__3.real = sgn * sumr.real;
                z__3.imag = sgn * sumr.imag; // , expr subst
                z__2.real = suml.real + z__3.real;
                z__2.imag = suml.imag + z__3.imag; // , expr subst
                z__1.real = c__[i__3].real - z__2.real;
                z__1.imag = c__[i__3].imag - z__2.imag; // , expr subst
                vec.real = z__1.real;
                vec.imag = z__1.imag; // , expr subst
                scaloc = 1.;
                d_cnjg(&z__2, &a[k + k * a_dim1]);
                i__3 = l + l * b_dim1;
                z__3.real = sgn * b[i__3].real;
                z__3.imag = sgn * b[i__3].imag; // , expr subst
                z__1.real = z__2.real + z__3.real;
                z__1.imag = z__2.imag + z__3.imag; // , expr subst
                a11.real = z__1.real;
                a11.imag = z__1.imag; // , expr subst
                da11 = (d__1 = a11.real, f2c_dabs(d__1)) + (d__2 = d_imag(&a11), f2c_dabs(d__2));
                if(da11 <= smin)
                {
                    a11.real = smin;
                    a11.imag = 0.; // , expr subst
                    da11 = smin;
                    *info = 1;
                }
                db = (d__1 = vec.real, f2c_dabs(d__1)) + (d__2 = d_imag(&vec), f2c_dabs(d__2));
                if(da11 < 1. && db > 1.)
                {
                    if(db > bignum * da11)
                    {
                        scaloc = 1. / db;
                    }
                }
                z__3.real = scaloc;
                z__3.imag = 0.; // , expr subst
                z__2.real = vec.real * z__3.real - vec.imag * z__3.imag;
                z__2.imag = vec.real * z__3.imag + vec.imag * z__3.real; // , expr subst
                zladiv_f2c_(&z__1, &z__2, &a11);
                x11.real = z__1.real;
                x11.imag = z__1.imag; // , expr subst
                if(scaloc != 1.)
                {
                    i__3 = *n;
                    for(j = 1; j <= i__3; ++j)
                    {
                        aocl_blas_zdscal(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                        /* L40: */
                    }
                    *scale *= scaloc;
                }
                i__3 = k + l * c_dim1;
                c__[i__3].real = x11.real;
                c__[i__3].imag = x11.imag; // , expr subst
                /* L50: */
            }
            /* L60: */
        }
    }
    else if(!notrna && !notrnb)
    {
        /* Solve A**H*X + ISGN*X*B**H = C. */
        /* The (K,L)th block of X is determined starting from */
        /* upper-right corner column by column by */
        /* A**H(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L) */
        /* Where */
        /* K-1 */
        /* R(K,L) = SUM [A**H(I,K)*X(I,L)] + */
        /* I=1 */
        /* N */
        /* ISGN*SUM [X(K,J)*B**H(L,J)]. */
        /* J=L+1 */
        for(l = *n; l >= 1; --l)
        {
            i__1 = *m;
            for(k = 1; k <= i__1; ++k)
            {
                i__2 = k - 1;
                aocl_lapack_zdotc_f2c(&z__1, &i__2, &a[k * a_dim1 + 1], &c__1, &c__[l * c_dim1 + 1], &c__1);
                suml.real = z__1.real;
                suml.imag = z__1.imag; // , expr subst
                i__2 = *n - l;
                /* Computing MIN */
                i__3 = l + 1;
                /* Computing MIN */
                i__4 = l + 1;
                aocl_lapack_zdotc_f2c(&z__1, &i__2, &c__[k + fla_min(i__3, *n) * c_dim1], ldc,
                           &b[l + fla_min(i__4, *n) * b_dim1], ldb);
                sumr.real = z__1.real;
                sumr.imag = z__1.imag; // , expr subst
                i__2 = k + l * c_dim1;
                d_cnjg(&z__4, &sumr);
                z__3.real = sgn * z__4.real;
                z__3.imag = sgn * z__4.imag; // , expr subst
                z__2.real = suml.real + z__3.real;
                z__2.imag = suml.imag + z__3.imag; // , expr subst
                z__1.real = c__[i__2].real - z__2.real;
                z__1.imag = c__[i__2].imag - z__2.imag; // , expr subst
                vec.real = z__1.real;
                vec.imag = z__1.imag; // , expr subst
                scaloc = 1.;
                i__2 = k + k * a_dim1;
                i__3 = l + l * b_dim1;
                z__3.real = sgn * b[i__3].real;
                z__3.imag = sgn * b[i__3].imag; // , expr subst
                z__2.real = a[i__2].real + z__3.real;
                z__2.imag = a[i__2].imag + z__3.imag; // , expr subst
                d_cnjg(&z__1, &z__2);
                a11.real = z__1.real;
                a11.imag = z__1.imag; // , expr subst
                da11 = (d__1 = a11.real, f2c_dabs(d__1)) + (d__2 = d_imag(&a11), f2c_dabs(d__2));
                if(da11 <= smin)
                {
                    a11.real = smin;
                    a11.imag = 0.; // , expr subst
                    da11 = smin;
                    *info = 1;
                }
                db = (d__1 = vec.real, f2c_dabs(d__1)) + (d__2 = d_imag(&vec), f2c_dabs(d__2));
                if(da11 < 1. && db > 1.)
                {
                    if(db > bignum * da11)
                    {
                        scaloc = 1. / db;
                    }
                }
                z__3.real = scaloc;
                z__3.imag = 0.; // , expr subst
                z__2.real = vec.real * z__3.real - vec.imag * z__3.imag;
                z__2.imag = vec.real * z__3.imag + vec.imag * z__3.real; // , expr subst
                zladiv_f2c_(&z__1, &z__2, &a11);
                x11.real = z__1.real;
                x11.imag = z__1.imag; // , expr subst
                if(scaloc != 1.)
                {
                    i__2 = *n;
                    for(j = 1; j <= i__2; ++j)
                    {
                        aocl_blas_zdscal(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                        /* L70: */
                    }
                    *scale *= scaloc;
                }
                i__2 = k + l * c_dim1;
                c__[i__2].real = x11.real;
                c__[i__2].imag = x11.imag; // , expr subst
                /* L80: */
            }
            /* L90: */
        }
    }
    else if(notrna && !notrnb)
    {
        /* Solve A*X + ISGN*X*B**H = C. */
        /* The (K,L)th block of X is determined starting from */
        /* bottom-left corner column by column by */
        /* A(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L) */
        /* Where */
        /* M N */
        /* R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B**H(L,J)] */
        /* I=K+1 J=L+1 */
        for(l = *n; l >= 1; --l)
        {
            for(k = *m; k >= 1; --k)
            {
                i__1 = *m - k;
                /* Computing MIN */
                i__2 = k + 1;
                /* Computing MIN */
                i__3 = k + 1;
                aocl_lapack_zdotu_f2c(&z__1, &i__1, &a[k + fla_min(i__2, *m) * a_dim1], lda,
                           &c__[fla_min(i__3, *m) + l * c_dim1], &c__1);
                suml.real = z__1.real;
                suml.imag = z__1.imag; // , expr subst
                i__1 = *n - l;
                /* Computing MIN */
                i__2 = l + 1;
                /* Computing MIN */
                i__3 = l + 1;
                aocl_lapack_zdotc_f2c(&z__1, &i__1, &c__[k + fla_min(i__2, *n) * c_dim1], ldc,
                           &b[l + fla_min(i__3, *n) * b_dim1], ldb);
                sumr.real = z__1.real;
                sumr.imag = z__1.imag; // , expr subst
                i__1 = k + l * c_dim1;
                d_cnjg(&z__4, &sumr);
                z__3.real = sgn * z__4.real;
                z__3.imag = sgn * z__4.imag; // , expr subst
                z__2.real = suml.real + z__3.real;
                z__2.imag = suml.imag + z__3.imag; // , expr subst
                z__1.real = c__[i__1].real - z__2.real;
                z__1.imag = c__[i__1].imag - z__2.imag; // , expr subst
                vec.real = z__1.real;
                vec.imag = z__1.imag; // , expr subst
                scaloc = 1.;
                i__1 = k + k * a_dim1;
                d_cnjg(&z__3, &b[l + l * b_dim1]);
                z__2.real = sgn * z__3.real;
                z__2.imag = sgn * z__3.imag; // , expr subst
                z__1.real = a[i__1].real + z__2.real;
                z__1.imag = a[i__1].imag + z__2.imag; // , expr subst
                a11.real = z__1.real;
                a11.imag = z__1.imag; // , expr subst
                da11 = (d__1 = a11.real, f2c_dabs(d__1)) + (d__2 = d_imag(&a11), f2c_dabs(d__2));
                if(da11 <= smin)
                {
                    a11.real = smin;
                    a11.imag = 0.; // , expr subst
                    da11 = smin;
                    *info = 1;
                }
                db = (d__1 = vec.real, f2c_dabs(d__1)) + (d__2 = d_imag(&vec), f2c_dabs(d__2));
                if(da11 < 1. && db > 1.)
                {
                    if(db > bignum * da11)
                    {
                        scaloc = 1. / db;
                    }
                }
                z__3.real = scaloc;
                z__3.imag = 0.; // , expr subst
                z__2.real = vec.real * z__3.real - vec.imag * z__3.imag;
                z__2.imag = vec.real * z__3.imag + vec.imag * z__3.real; // , expr subst
                zladiv_f2c_(&z__1, &z__2, &a11);
                x11.real = z__1.real;
                x11.imag = z__1.imag; // , expr subst
                if(scaloc != 1.)
                {
                    i__1 = *n;
                    for(j = 1; j <= i__1; ++j)
                    {
                        aocl_blas_zdscal(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
                        /* L100: */
                    }
                    *scale *= scaloc;
                }
                i__1 = k + l * c_dim1;
                c__[i__1].real = x11.real;
                c__[i__1].imag = x11.imag; // , expr subst
                /* L110: */
            }
            /* L120: */
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZTRSYL */
}
/* ztrsyl_ */
