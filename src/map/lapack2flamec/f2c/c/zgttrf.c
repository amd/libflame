/* ../netlib/zgttrf.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZGTTRF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGTTRF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgttrf.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgttrf.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgttrf.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGTTRF( N, DL, D, DU, DU2, IPIV, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 D( * ), DL( * ), DU( * ), DU2( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGTTRF computes an LU factorization of a scomplex tridiagonal matrix A */
/* > using elimination with partial pivoting and row interchanges. */
/* > */
/* > The factorization has the form */
/* > A = L * U */
/* > where L is a product of permutation and unit lower bidiagonal */
/* > matrices and U is upper triangular with nonzeros in only the main */
/* > diagonal and first two superdiagonals. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DL */
/* > \verbatim */
/* > DL is COMPLEX*16 array, dimension (N-1) */
/* > On entry, DL must contain the (n-1) sub-diagonal elements of */
/* > A. */
/* > */
/* > On exit, DL is overwritten by the (n-1) multipliers that */
/* > define the matrix L from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is COMPLEX*16 array, dimension (N) */
/* > On entry, D must contain the diagonal elements of A. */
/* > */
/* > On exit, D is overwritten by the n diagonal elements of the */
/* > upper triangular matrix U from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU */
/* > \verbatim */
/* > DU is COMPLEX*16 array, dimension (N-1) */
/* > On entry, DU must contain the (n-1) super-diagonal elements */
/* > of A. */
/* > */
/* > On exit, DU is overwritten by the (n-1) elements of the first */
/* > super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[out] DU2 */
/* > \verbatim */
/* > DU2 is COMPLEX*16 array, dimension (N-2) */
/* > On exit, DU2 is overwritten by the (n-2) elements of the */
/* > second super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > The pivot indices;
for 1 <= i <= n, row i of the matrix was */
/* > interchanged with row IPIV(i). IPIV(i) will always be either */
/* > i or i+1;
IPIV(i) = i indicates a row interchange was not */
/* > required. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -k, the k-th argument had an illegal value */
/* > > 0: if INFO = k, U(k,k) is exactly zero. The factorization */
/* > has been completed, but the factor U is exactly */
/* > singular, and division by zero will occur if it is used */
/* > to solve a system of equations. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16GTcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zgttrf_(aocl_int_t *n, dcomplex *dl, dcomplex *d__, dcomplex *du,
             dcomplex *du2, aocl_int_t *ipiv, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zgttrf(n, dl, d__, du, du2, ipiv, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zgttrf(&n_64, dl, d__, du, du2, ipiv, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_zgttrf(aocl_int64_t *n, dcomplex *dl, dcomplex *d__, dcomplex *du,
                        dcomplex *du2, aocl_int_t *ipiv, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgttrf inputs: n %" FLA_IS "", *n);

    /* System generated locals */
    aocl_int64_t i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;
    dcomplex z__1, z__2;
    /* Builtin functions */
    double d_imag(dcomplex *);
    void z_div(dcomplex *, dcomplex *, dcomplex *);
    /* Local variables */
    aocl_int64_t i__;
    dcomplex fact, temp;
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --ipiv;
    --du2;
    --du;
    --d__;
    --dl;
    /* Function Body */
    *info = 0;
    if(*n < 0)
    {
        *info = -1;
        i__1 = -(*info);
        aocl_blas_xerbla("ZGTTRF", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Initialize IPIV(i) = i and DU2(i) = 0 */
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        ipiv[i__] = (aocl_int_t)(i__);
        /* L10: */
    }
    i__1 = *n - 2;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = i__;
        du2[i__2].real = 0.;
        du2[i__2].imag = 0.; // , expr subst
        /* L20: */
    }
    i__1 = *n - 2;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = i__;
        i__3 = i__;
        if((d__1 = d__[i__2].real, f2c_dabs(d__1)) + (d__2 = d_imag(&d__[i__]), f2c_dabs(d__2))
           >= (d__3 = dl[i__3].real, f2c_dabs(d__3)) + (d__4 = d_imag(&dl[i__]), f2c_dabs(d__4)))
        {
            /* No row interchange required, eliminate DL(I) */
            i__2 = i__;
            if((d__1 = d__[i__2].real, f2c_dabs(d__1)) + (d__2 = d_imag(&d__[i__]), f2c_dabs(d__2))
               != 0.)
            {
                z_div(&z__1, &dl[i__], &d__[i__]);
                fact.real = z__1.real;
                fact.imag = z__1.imag; // , expr subst
                i__2 = i__;
                dl[i__2].real = fact.real;
                dl[i__2].imag = fact.imag; // , expr subst
                i__2 = i__ + 1;
                i__3 = i__ + 1;
                i__4 = i__;
                z__2.real = fact.real * du[i__4].real - fact.imag * du[i__4].imag;
                z__2.imag = fact.real * du[i__4].imag + fact.imag * du[i__4].real; // , expr subst
                z__1.real = d__[i__3].real - z__2.real;
                z__1.imag = d__[i__3].imag - z__2.imag; // , expr subst
                d__[i__2].real = z__1.real;
                d__[i__2].imag = z__1.imag; // , expr subst
            }
        }
        else
        {
            /* Interchange rows I and I+1, eliminate DL(I) */
            z_div(&z__1, &d__[i__], &dl[i__]);
            fact.real = z__1.real;
            fact.imag = z__1.imag; // , expr subst
            i__2 = i__;
            i__3 = i__;
            d__[i__2].real = dl[i__3].real;
            d__[i__2].imag = dl[i__3].imag; // , expr subst
            i__2 = i__;
            dl[i__2].real = fact.real;
            dl[i__2].imag = fact.imag; // , expr subst
            i__2 = i__;
            temp.real = du[i__2].real;
            temp.imag = du[i__2].imag; // , expr subst
            i__2 = i__;
            i__3 = i__ + 1;
            du[i__2].real = d__[i__3].real;
            du[i__2].imag = d__[i__3].imag; // , expr subst
            i__2 = i__ + 1;
            i__3 = i__ + 1;
            z__2.real = fact.real * d__[i__3].real - fact.imag * d__[i__3].imag;
            z__2.imag = fact.real * d__[i__3].imag + fact.imag * d__[i__3].real; // , expr subst
            z__1.real = temp.real - z__2.real;
            z__1.imag = temp.imag - z__2.imag; // , expr subst
            d__[i__2].real = z__1.real;
            d__[i__2].imag = z__1.imag; // , expr subst
            i__2 = i__;
            i__3 = i__ + 1;
            du2[i__2].real = du[i__3].real;
            du2[i__2].imag = du[i__3].imag; // , expr subst
            i__2 = i__ + 1;
            z__2.real = -fact.real;
            z__2.imag = -fact.imag; // , expr subst
            i__3 = i__ + 1;
            z__1.real = z__2.real * du[i__3].real - z__2.imag * du[i__3].imag;
            z__1.imag = z__2.real * du[i__3].imag + z__2.imag * du[i__3].real; // , expr subst
            du[i__2].real = z__1.real;
            du[i__2].imag = z__1.imag; // , expr subst
            ipiv[i__] = (aocl_int_t)(i__ + 1);
        }
        /* L30: */
    }
    if(*n > 1)
    {
        i__ = *n - 1;
        i__1 = i__;
        i__2 = i__;
        if((d__1 = d__[i__1].real, f2c_dabs(d__1)) + (d__2 = d_imag(&d__[i__]), f2c_dabs(d__2))
           >= (d__3 = dl[i__2].real, f2c_dabs(d__3)) + (d__4 = d_imag(&dl[i__]), f2c_dabs(d__4)))
        {
            i__1 = i__;
            if((d__1 = d__[i__1].real, f2c_dabs(d__1)) + (d__2 = d_imag(&d__[i__]), f2c_dabs(d__2))
               != 0.)
            {
                z_div(&z__1, &dl[i__], &d__[i__]);
                fact.real = z__1.real;
                fact.imag = z__1.imag; // , expr subst
                i__1 = i__;
                dl[i__1].real = fact.real;
                dl[i__1].imag = fact.imag; // , expr subst
                i__1 = i__ + 1;
                i__2 = i__ + 1;
                i__3 = i__;
                z__2.real = fact.real * du[i__3].real - fact.imag * du[i__3].imag;
                z__2.imag = fact.real * du[i__3].imag + fact.imag * du[i__3].real; // , expr subst
                z__1.real = d__[i__2].real - z__2.real;
                z__1.imag = d__[i__2].imag - z__2.imag; // , expr subst
                d__[i__1].real = z__1.real;
                d__[i__1].imag = z__1.imag; // , expr subst
            }
        }
        else
        {
            z_div(&z__1, &d__[i__], &dl[i__]);
            fact.real = z__1.real;
            fact.imag = z__1.imag; // , expr subst
            i__1 = i__;
            i__2 = i__;
            d__[i__1].real = dl[i__2].real;
            d__[i__1].imag = dl[i__2].imag; // , expr subst
            i__1 = i__;
            dl[i__1].real = fact.real;
            dl[i__1].imag = fact.imag; // , expr subst
            i__1 = i__;
            temp.real = du[i__1].real;
            temp.imag = du[i__1].imag; // , expr subst
            i__1 = i__;
            i__2 = i__ + 1;
            du[i__1].real = d__[i__2].real;
            du[i__1].imag = d__[i__2].imag; // , expr subst
            i__1 = i__ + 1;
            i__2 = i__ + 1;
            z__2.real = fact.real * d__[i__2].real - fact.imag * d__[i__2].imag;
            z__2.imag = fact.real * d__[i__2].imag + fact.imag * d__[i__2].real; // , expr subst
            z__1.real = temp.real - z__2.real;
            z__1.imag = temp.imag - z__2.imag; // , expr subst
            d__[i__1].real = z__1.real;
            d__[i__1].imag = z__1.imag; // , expr subst
            ipiv[i__] = (aocl_int_t)(i__ + 1);
        }
    }
    /* Check for a zero on the diagonal of U. */
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = i__;
        if((d__1 = d__[i__2].real, f2c_dabs(d__1)) + (d__2 = d_imag(&d__[i__]), f2c_dabs(d__2)) == 0.)
        {
            *info = i__;
            goto L50;
        }
        /* L40: */
    }
L50:
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZGTTRF */
}
/* zgttrf_ */
