/* ../netlib/zgtsv.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief <b> ZGTSV computes the solution to system of linear equations A * X = B for GT matrices <b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGTSV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgtsv.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgtsv.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgtsv.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGTSV( N, NRHS, DL, D, DU, B, LDB, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 B( LDB, * ), D( * ), DL( * ), DU( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGTSV solves the equation */
/* > */
/* > A*X = B, */
/* > */
/* > where A is an N-by-N tridiagonal matrix, by Gaussian elimination with */
/* > partial pivoting. */
/* > */
/* > Note that the equation A**T *X = B may be solved by interchanging the */
/* > order of the arguments DU and DL. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DL */
/* > \verbatim */
/* > DL is COMPLEX*16 array, dimension (N-1) */
/* > On entry, DL must contain the (n-1) subdiagonal elements of */
/* > A. */
/* > On exit, DL is overwritten by the (n-2) elements of the */
/* > second superdiagonal of the upper triangular matrix U from */
/* > the LU factorization of A, in DL(1), ..., DL(n-2). */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is COMPLEX*16 array, dimension (N) */
/* > On entry, D must contain the diagonal elements of A. */
/* > On exit, D is overwritten by the n diagonal elements of U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU */
/* > \verbatim */
/* > DU is COMPLEX*16 array, dimension (N-1) */
/* > On entry, DU must contain the (n-1) superdiagonal elements */
/* > of A. */
/* > On exit, DU is overwritten by the (n-1) elements of the first */
/* > superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* > On entry, the N-by-NRHS right hand side matrix B. */
/* > On exit, if INFO = 0, the N-by-NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, U(i,i) is exactly zero, and the solution */
/* > has not been computed. The factorization has not been */
/* > completed unless i = N. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16GTsolve */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zgtsv_(aocl_int_t *n, aocl_int_t *nrhs, dcomplex *dl, dcomplex *d__,
            dcomplex *du, dcomplex *b, aocl_int_t *ldb, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zgtsv(n, nrhs, dl, d__, du, b, ldb, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zgtsv(&n_64, &nrhs_64, dl, d__, du, b, &ldb_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_zgtsv(aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *dl, dcomplex *d__,
                       dcomplex *du, dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgtsv inputs: n %" FLA_IS ", nrhs %" FLA_IS ", ldb %" FLA_IS "", *n, *nrhs,
                      *ldb);

    /* System generated locals */
    aocl_int64_t b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublereal d__1, d__2, d__3, d__4;
    dcomplex z__1, z__2, z__3, z__4, z__5;
    /* Builtin functions */
    double d_imag(dcomplex *);
    void z_div(dcomplex *, dcomplex *, dcomplex *);
    /* Local variables */
    aocl_int64_t j, k;
    dcomplex temp, mult;
    /* -- LAPACK driver routine (version 3.4.2) -- */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --dl;
    --d__;
    --du;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    *info = 0;
    if(*n < 0)
    {
        *info = -1;
    }
    else if(*nrhs < 0)
    {
        *info = -2;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -7;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("ZGTSV ", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    i__1 = *n - 1;
    for(k = 1; k <= i__1; ++k)
    {
        i__2 = k;
        if(dl[i__2].real == 0. && dl[i__2].imag == 0.)
        {
            /* Subdiagonal is zero, no elimination is required. */
            i__2 = k;
            if(d__[i__2].real == 0. && d__[i__2].imag == 0.)
            {
                /* Diagonal is zero: set INFO = K and return;
                a unique */
                /* solution can not be found. */
                *info = k;
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
        }
        else /* if(complicated condition) */
        {
            i__2 = k;
            i__3 = k;
            if((d__1 = d__[i__2].real, f2c_dabs(d__1)) + (d__2 = d_imag(&d__[k]), f2c_dabs(d__2))
               >= (d__3 = dl[i__3].real, f2c_dabs(d__3)) + (d__4 = d_imag(&dl[k]), f2c_dabs(d__4)))
            {
                /* No row interchange required */
                z_div(&z__1, &dl[k], &d__[k]);
                mult.real = z__1.real;
                mult.imag = z__1.imag; // , expr subst
                i__2 = k + 1;
                i__3 = k + 1;
                i__4 = k;
                z__2.real = mult.real * du[i__4].real - mult.imag * du[i__4].imag;
                z__2.imag = mult.real * du[i__4].imag + mult.imag * du[i__4].real; // , expr subst
                z__1.real = d__[i__3].real - z__2.real;
                z__1.imag = d__[i__3].imag - z__2.imag; // , expr subst
                d__[i__2].real = z__1.real;
                d__[i__2].imag = z__1.imag; // , expr subst
                i__2 = *nrhs;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = k + 1 + j * b_dim1;
                    i__4 = k + 1 + j * b_dim1;
                    i__5 = k + j * b_dim1;
                    z__2.real = mult.real * b[i__5].real - mult.imag * b[i__5].imag;
                    z__2.imag = mult.real * b[i__5].imag + mult.imag * b[i__5].real; // , expr subst
                    z__1.real = b[i__4].real - z__2.real;
                    z__1.imag = b[i__4].imag - z__2.imag; // , expr subst
                    b[i__3].real = z__1.real;
                    b[i__3].imag = z__1.imag; // , expr subst
                    /* L10: */
                }
                if(k < *n - 1)
                {
                    i__2 = k;
                    dl[i__2].real = 0.;
                    dl[i__2].imag = 0.; // , expr subst
                }
            }
            else
            {
                /* Interchange rows K and K+1 */
                z_div(&z__1, &d__[k], &dl[k]);
                mult.real = z__1.real;
                mult.imag = z__1.imag; // , expr subst
                i__2 = k;
                i__3 = k;
                d__[i__2].real = dl[i__3].real;
                d__[i__2].imag = dl[i__3].imag; // , expr subst
                i__2 = k + 1;
                temp.real = d__[i__2].real;
                temp.imag = d__[i__2].imag; // , expr subst
                i__2 = k + 1;
                i__3 = k;
                z__2.real = mult.real * temp.real - mult.imag * temp.imag;
                z__2.imag = mult.real * temp.imag + mult.imag * temp.real; // , expr subst
                z__1.real = du[i__3].real - z__2.real;
                z__1.imag = du[i__3].imag - z__2.imag; // , expr subst
                d__[i__2].real = z__1.real;
                d__[i__2].imag = z__1.imag; // , expr subst
                if(k < *n - 1)
                {
                    i__2 = k;
                    i__3 = k + 1;
                    dl[i__2].real = du[i__3].real;
                    dl[i__2].imag = du[i__3].imag; // , expr subst
                    i__2 = k + 1;
                    z__2.real = -mult.real;
                    z__2.imag = -mult.imag; // , expr subst
                    i__3 = k;
                    z__1.real = z__2.real * dl[i__3].real - z__2.imag * dl[i__3].imag;
                    z__1.imag = z__2.real * dl[i__3].imag + z__2.imag * dl[i__3].real; // , expr subst
                    du[i__2].real = z__1.real;
                    du[i__2].imag = z__1.imag; // , expr subst
                }
                i__2 = k;
                du[i__2].real = temp.real;
                du[i__2].imag = temp.imag; // , expr subst
                i__2 = *nrhs;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = k + j * b_dim1;
                    temp.real = b[i__3].real;
                    temp.imag = b[i__3].imag; // , expr subst
                    i__3 = k + j * b_dim1;
                    i__4 = k + 1 + j * b_dim1;
                    b[i__3].real = b[i__4].real;
                    b[i__3].imag = b[i__4].imag; // , expr subst
                    i__3 = k + 1 + j * b_dim1;
                    i__4 = k + 1 + j * b_dim1;
                    z__2.real = mult.real * b[i__4].real - mult.imag * b[i__4].imag;
                    z__2.imag = mult.real * b[i__4].imag + mult.imag * b[i__4].real; // , expr subst
                    z__1.real = temp.real - z__2.real;
                    z__1.imag = temp.imag - z__2.imag; // , expr subst
                    b[i__3].real = z__1.real;
                    b[i__3].imag = z__1.imag; // , expr subst
                    /* L20: */
                }
            }
        }
        /* L30: */
    }
    i__1 = *n;
    if(d__[i__1].real == 0. && d__[i__1].imag == 0.)
    {
        *info = *n;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Back solve with the matrix U from the factorization. */
    i__1 = *nrhs;
    for(j = 1; j <= i__1; ++j)
    {
        i__2 = *n + j * b_dim1;
        z_div(&z__1, &b[*n + j * b_dim1], &d__[*n]);
        b[i__2].real = z__1.real;
        b[i__2].imag = z__1.imag; // , expr subst
        if(*n > 1)
        {
            i__2 = *n - 1 + j * b_dim1;
            i__3 = *n - 1 + j * b_dim1;
            i__4 = *n - 1;
            i__5 = *n + j * b_dim1;
            z__3.real = du[i__4].real * b[i__5].real - du[i__4].imag * b[i__5].imag;
            z__3.imag = du[i__4].real * b[i__5].imag + du[i__4].imag * b[i__5].real; // , expr subst
            z__2.real = b[i__3].real - z__3.real;
            z__2.imag = b[i__3].imag - z__3.imag; // , expr subst
            z_div(&z__1, &z__2, &d__[*n - 1]);
            b[i__2].real = z__1.real;
            b[i__2].imag = z__1.imag; // , expr subst
        }
        for(k = *n - 2; k >= 1; --k)
        {
            i__2 = k + j * b_dim1;
            i__3 = k + j * b_dim1;
            i__4 = k;
            i__5 = k + 1 + j * b_dim1;
            z__4.real = du[i__4].real * b[i__5].real - du[i__4].imag * b[i__5].imag;
            z__4.imag = du[i__4].real * b[i__5].imag + du[i__4].imag * b[i__5].real; // , expr subst
            z__3.real = b[i__3].real - z__4.real;
            z__3.imag = b[i__3].imag - z__4.imag; // , expr subst
            i__6 = k;
            i__7 = k + 2 + j * b_dim1;
            z__5.real = dl[i__6].real * b[i__7].real - dl[i__6].imag * b[i__7].imag;
            z__5.imag = dl[i__6].real * b[i__7].imag + dl[i__6].imag * b[i__7].real; // , expr subst
            z__2.real = z__3.real - z__5.real;
            z__2.imag = z__3.imag - z__5.imag; // , expr subst
            z_div(&z__1, &z__2, &d__[k]);
            b[i__2].real = z__1.real;
            b[i__2].imag = z__1.imag; // , expr subst
            /* L40: */
        }
        /* L50: */
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZGTSV */
}
/* zgtsv_ */
