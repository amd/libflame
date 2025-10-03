/* ./zgesc2.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
static dcomplex c_b13 = {1., 0.};
static aocl_int64_t c_n1 = -1;
/* > \brief \b ZGESC2 solves a system of linear equations using the LU factorization with complete
 * pivoting co mputed by sgetc2. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGESC2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgesc2.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgesc2.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgesc2.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE ) */
/* .. Scalar Arguments .. */
/* INTEGER LDA, N */
/* DOUBLE PRECISION SCALE */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ), JPIV( * ) */
/* COMPLEX*16 A( LDA, * ), RHS( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGESC2 solves a system of linear equations */
/* > */
/* > A * X = scale* RHS */
/* > */
/* > with a general N-by-N matrix A using the LU factorization with */
/* > complete pivoting computed by ZGETC2. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA, N) */
/* > On entry, the LU part of the factorization of the n-by-n */
/* > matrix A computed by ZGETC2: A = P * L * U * Q */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] RHS */
/* > \verbatim */
/* > RHS is COMPLEX*16 array, dimension N. */
/* > On entry, the right hand side vector b. */
/* > On exit, the solution vector X. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N). */
/* > The pivot indices;
for 1 <= i <= N, row i of the */
/* > matrix has been interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] JPIV */
/* > \verbatim */
/* > JPIV is INTEGER array, dimension (N). */
/* > The pivot indices;
for 1 <= j <= N, column j of the */
/* > matrix has been interchanged with column JPIV(j). */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* > SCALE is DOUBLE PRECISION */
/* > On exit, SCALE contains the scale factor. SCALE is chosen */
/* > 0 <= SCALE <= 1 to prevent overflow in the solution. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup gesc2 */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* > Umea University, S-901 87 Umea, Sweden. */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zgesc2_(aocl_int_t *n, dcomplex *a, aocl_int_t *lda, dcomplex *rhs, aocl_int_t *ipiv,
             aocl_int_t *jpiv, doublereal *scale)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zgesc2(n, a, lda, rhs, ipiv, jpiv, scale);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;

    aocl_lapack_zgesc2(&n_64, a, &lda_64, rhs, ipiv, jpiv, scale);
#endif
}

void aocl_lapack_zgesc2(aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *rhs,
                        aocl_int_t *ipiv, aocl_int_t *jpiv, doublereal *scale)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgesc2 inputs: n %" FLA_IS ", lda %" FLA_IS "", *n, *lda);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    dcomplex z__1, z__2, z__3;
    /* Builtin functions */
    double z_abs(dcomplex *);
    void z_div(dcomplex *, dcomplex *, dcomplex *);
    /* Local variables */
    aocl_int64_t i__, j;
    doublereal eps;
    dcomplex temp;
    extern doublereal dlamch_(char *);
    doublereal smlnum;
    /* -- LAPACK auxiliary routine -- */
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Set constant to control overflow */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --rhs;
    --ipiv;
    --jpiv;
    /* Function Body */
    eps = dlamch_("P");
    smlnum = dlamch_("S") / eps;
    /* Apply permutations IPIV to RHS */
    i__1 = *n - 1;
    aocl_lapack_zlaswp(&c__1, &rhs[1], lda, &c__1, &i__1, &ipiv[1], &c__1);
    /* Solve for L part */
    i__1 = *n - 1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = *n;
        for(j = i__ + 1; j <= i__2; ++j)
        {
            i__3 = j;
            i__4 = j;
            i__5 = j + i__ * a_dim1;
            i__6 = i__;
            z__2.real = a[i__5].real * rhs[i__6].real - a[i__5].imag * rhs[i__6].imag;
            z__2.imag = a[i__5].real * rhs[i__6].imag + a[i__5].imag * rhs[i__6].real; // , expr subst
            z__1.real = rhs[i__4].real - z__2.real;
            z__1.imag = rhs[i__4].imag - z__2.imag; // , expr subst
            rhs[i__3].real = z__1.real;
            rhs[i__3].imag = z__1.imag; // , expr subst
            /* L10: */
        }
        /* L20: */
    }
    /* Solve for U part */
    *scale = 1.;
    /* Check for scaling */
    i__ = aocl_blas_izamax(n, &rhs[1], &c__1);
    if(smlnum * 2. * z_abs(&rhs[i__]) > z_abs(&a[*n + *n * a_dim1]))
    {
        d__1 = z_abs(&rhs[i__]);
        z__1.real = .5 / d__1;
        z__1.imag = 0. / d__1; // , expr subst
        temp.real = z__1.real;
        temp.imag = z__1.imag; // , expr subst
        aocl_blas_zscal(n, &temp, &rhs[1], &c__1);
        *scale *= temp.real;
    }
    for(i__ = *n; i__ >= 1; --i__)
    {
        z_div(&z__1, &c_b13, &a[i__ + i__ * a_dim1]);
        temp.real = z__1.real;
        temp.imag = z__1.imag; // , expr subst
        i__1 = i__;
        i__2 = i__;
        z__1.real = rhs[i__2].real * temp.real - rhs[i__2].imag * temp.imag;
        z__1.imag = rhs[i__2].real * temp.imag + rhs[i__2].imag * temp.real; // , expr subst
        rhs[i__1].real = z__1.real;
        rhs[i__1].imag = z__1.imag; // , expr subst
        i__1 = *n;
        for(j = i__ + 1; j <= i__1; ++j)
        {
            i__2 = i__;
            i__3 = i__;
            i__4 = j;
            i__5 = i__ + j * a_dim1;
            z__3.real = a[i__5].real * temp.real - a[i__5].imag * temp.imag;
            z__3.imag = a[i__5].real * temp.imag + a[i__5].imag * temp.real; // , expr subst
            z__2.real = rhs[i__4].real * z__3.real - rhs[i__4].imag * z__3.imag;
            z__2.imag = rhs[i__4].real * z__3.imag + rhs[i__4].imag * z__3.real; // , expr subst
            z__1.real = rhs[i__3].real - z__2.real;
            z__1.imag = rhs[i__3].imag - z__2.imag; // , expr subst
            rhs[i__2].real = z__1.real;
            rhs[i__2].imag = z__1.imag; // , expr subst
            /* L30: */
        }
        /* L40: */
    }
    /* Apply permutations JPIV to the solution (RHS) */
    i__1 = *n - 1;
    aocl_lapack_zlaswp(&c__1, &rhs[1], lda, &c__1, &i__1, &jpiv[1], &c_n1);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZGESC2 */
}
/* zgesc2_ */
