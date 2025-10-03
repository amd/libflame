/* ../netlib/zgtts2.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZGTTS2 solves a system of linear equations with a tridiagonal matrix using the LU factorization computed by sgttrf. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGTTS2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgtts2.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgtts2.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgtts2.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB ) */
/* .. Scalar Arguments .. */
/* INTEGER ITRANS, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGTTS2 solves one of the systems of equations */
/* > A * X = B, A**T * X = B, or A**H * X = B, */
/* > with a tridiagonal matrix A using the LU factorization computed */
/* > by ZGTTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ITRANS */
/* > \verbatim */
/* > ITRANS is INTEGER */
/* > Specifies the form of the system of equations. */
/* > = 0: A * X = B (No transpose) */
/* > = 1: A**T * X = B (Transpose) */
/* > = 2: A**H * X = B (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* > DL is COMPLEX*16 array, dimension (N-1) */
/* > The (n-1) multipliers that define the matrix L from the */
/* > LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is COMPLEX*16 array, dimension (N) */
/* > The n diagonal elements of the upper triangular matrix U from */
/* > the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* > DU is COMPLEX*16 array, dimension (N-1) */
/* > The (n-1) elements of the first super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* > DU2 is COMPLEX*16 array, dimension (N-2) */
/* > The (n-2) elements of the second super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
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
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* > On entry, the matrix of right hand side vectors B. */
/* > On exit, B is overwritten by the solution vectors X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
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
void zgtts2_(aocl_int_t *itrans, aocl_int_t *n, aocl_int_t *nrhs, dcomplex *dl,
             dcomplex *d__, dcomplex *du, dcomplex *du2, aocl_int_t *ipiv,
             dcomplex *b, aocl_int_t *ldb)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zgtts2(itrans, n, nrhs, dl, d__, du, du2, ipiv, b, ldb);
#else
    aocl_int64_t itrans_64 = *itrans;
    aocl_int64_t n_64 = *n;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t ldb_64 = *ldb;

    aocl_lapack_zgtts2(&itrans_64, &n_64, &nrhs_64, dl, d__, du, du2, ipiv, b, &ldb_64);
#endif
}

void aocl_lapack_zgtts2(aocl_int64_t *itrans, aocl_int64_t *n, aocl_int64_t *nrhs,
                        dcomplex *dl, dcomplex *d__, dcomplex *du,
                        dcomplex *du2, aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgtts2 inputs: itrans %" FLA_IS ", n %" FLA_IS ", nrhs %" FLA_IS
                      ", ldb %" FLA_IS "",
                      *itrans, *n, *nrhs, *ldb);

    /* System generated locals */
    aocl_int64_t b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    dcomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8;
    /* Builtin functions */
    void z_div(dcomplex *, dcomplex *, dcomplex *),
        d_cnjg(dcomplex *, dcomplex *);
    /* Local variables */
    aocl_int64_t i__, j;
    dcomplex temp;
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    --dl;
    --d__;
    --du;
    --du2;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    if(*n == 0 || *nrhs == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*itrans == 0)
    {
        /* Solve A*X = B using the LU factorization of A, */
        /* overwriting each right hand side vector with its solution. */
        if(*nrhs <= 1)
        {
            j = 1;
        L10: /* Solve L*x = b. */
            i__1 = *n - 1;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                if(ipiv[i__] == i__)
                {
                    i__2 = i__ + 1 + j * b_dim1;
                    i__3 = i__ + 1 + j * b_dim1;
                    i__4 = i__;
                    i__5 = i__ + j * b_dim1;
                    z__2.real = dl[i__4].real * b[i__5].real - dl[i__4].imag * b[i__5].imag;
                    z__2.imag = dl[i__4].real * b[i__5].imag + dl[i__4].imag * b[i__5].real; // , expr subst
                    z__1.real = b[i__3].real - z__2.real;
                    z__1.imag = b[i__3].imag - z__2.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                }
                else
                {
                    i__2 = i__ + j * b_dim1;
                    temp.real = b[i__2].real;
                    temp.imag = b[i__2].imag; // , expr subst
                    i__2 = i__ + j * b_dim1;
                    i__3 = i__ + 1 + j * b_dim1;
                    b[i__2].real = b[i__3].real;
                    b[i__2].imag = b[i__3].imag; // , expr subst
                    i__2 = i__ + 1 + j * b_dim1;
                    i__3 = i__;
                    i__4 = i__ + j * b_dim1;
                    z__2.real = dl[i__3].real * b[i__4].real - dl[i__3].imag * b[i__4].imag;
                    z__2.imag = dl[i__3].real * b[i__4].imag + dl[i__3].imag * b[i__4].real; // , expr subst
                    z__1.real = temp.real - z__2.real;
                    z__1.imag = temp.imag - z__2.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                }
                /* L20: */
            }
            /* Solve U*x = b. */
            i__1 = *n + j * b_dim1;
            z_div(&z__1, &b[*n + j * b_dim1], &d__[*n]);
            b[i__1].real = z__1.real;
            b[i__1].imag = z__1.imag; // , expr subst
            if(*n > 1)
            {
                i__1 = *n - 1 + j * b_dim1;
                i__2 = *n - 1 + j * b_dim1;
                i__3 = *n - 1;
                i__4 = *n + j * b_dim1;
                z__3.real = du[i__3].real * b[i__4].real - du[i__3].imag * b[i__4].imag;
                z__3.imag = du[i__3].real * b[i__4].imag + du[i__3].imag * b[i__4].real; // , expr subst
                z__2.real = b[i__2].real - z__3.real;
                z__2.imag = b[i__2].imag - z__3.imag; // , expr subst
                z_div(&z__1, &z__2, &d__[*n - 1]);
                b[i__1].real = z__1.real;
                b[i__1].imag = z__1.imag; // , expr subst
            }
            for(i__ = *n - 2; i__ >= 1; --i__)
            {
                i__1 = i__ + j * b_dim1;
                i__2 = i__ + j * b_dim1;
                i__3 = i__;
                i__4 = i__ + 1 + j * b_dim1;
                z__4.real = du[i__3].real * b[i__4].real - du[i__3].imag * b[i__4].imag;
                z__4.imag = du[i__3].real * b[i__4].imag + du[i__3].imag * b[i__4].real; // , expr subst
                z__3.real = b[i__2].real - z__4.real;
                z__3.imag = b[i__2].imag - z__4.imag; // , expr subst
                i__5 = i__;
                i__6 = i__ + 2 + j * b_dim1;
                z__5.real = du2[i__5].real * b[i__6].real - du2[i__5].imag * b[i__6].imag;
                z__5.imag = du2[i__5].real * b[i__6].imag + du2[i__5].imag * b[i__6].real; // , expr subst
                z__2.real = z__3.real - z__5.real;
                z__2.imag = z__3.imag - z__5.imag; // , expr subst
                z_div(&z__1, &z__2, &d__[i__]);
                b[i__1].real = z__1.real;
                b[i__1].imag = z__1.imag; // , expr subst
                /* L30: */
            }
            if(j < *nrhs)
            {
                ++j;
                goto L10;
            }
        }
        else
        {
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                /* Solve L*x = b. */
                i__2 = *n - 1;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    if(ipiv[i__] == i__)
                    {
                        i__3 = i__ + 1 + j * b_dim1;
                        i__4 = i__ + 1 + j * b_dim1;
                        i__5 = i__;
                        i__6 = i__ + j * b_dim1;
                        z__2.real = dl[i__5].real * b[i__6].real - dl[i__5].imag * b[i__6].imag;
                        z__2.imag = dl[i__5].real * b[i__6].imag + dl[i__5].imag * b[i__6].real; // , expr subst
                        z__1.real = b[i__4].real - z__2.real;
                        z__1.imag = b[i__4].imag - z__2.imag; // , expr subst
                        b[i__3].real = z__1.real;
                        b[i__3].imag = z__1.imag; // , expr subst
                    }
                    else
                    {
                        i__3 = i__ + j * b_dim1;
                        temp.real = b[i__3].real;
                        temp.imag = b[i__3].imag; // , expr subst
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + 1 + j * b_dim1;
                        b[i__3].real = b[i__4].real;
                        b[i__3].imag = b[i__4].imag; // , expr subst
                        i__3 = i__ + 1 + j * b_dim1;
                        i__4 = i__;
                        i__5 = i__ + j * b_dim1;
                        z__2.real = dl[i__4].real * b[i__5].real - dl[i__4].imag * b[i__5].imag;
                        z__2.imag = dl[i__4].real * b[i__5].imag + dl[i__4].imag * b[i__5].real; // , expr subst
                        z__1.real = temp.real - z__2.real;
                        z__1.imag = temp.imag - z__2.imag; // , expr subst
                        b[i__3].real = z__1.real;
                        b[i__3].imag = z__1.imag; // , expr subst
                    }
                    /* L40: */
                }
                /* Solve U*x = b. */
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
                for(i__ = *n - 2; i__ >= 1; --i__)
                {
                    i__2 = i__ + j * b_dim1;
                    i__3 = i__ + j * b_dim1;
                    i__4 = i__;
                    i__5 = i__ + 1 + j * b_dim1;
                    z__4.real = du[i__4].real * b[i__5].real - du[i__4].imag * b[i__5].imag;
                    z__4.imag = du[i__4].real * b[i__5].imag + du[i__4].imag * b[i__5].real; // , expr subst
                    z__3.real = b[i__3].real - z__4.real;
                    z__3.imag = b[i__3].imag - z__4.imag; // , expr subst
                    i__6 = i__;
                    i__7 = i__ + 2 + j * b_dim1;
                    z__5.real = du2[i__6].real * b[i__7].real - du2[i__6].imag * b[i__7].imag;
                    z__5.imag = du2[i__6].real * b[i__7].imag + du2[i__6].imag * b[i__7].real; // , expr subst
                    z__2.real = z__3.real - z__5.real;
                    z__2.imag = z__3.imag - z__5.imag; // , expr subst
                    z_div(&z__1, &z__2, &d__[i__]);
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                    /* L50: */
                }
                /* L60: */
            }
        }
    }
    else if(*itrans == 1)
    {
        /* Solve A**T * X = B. */
        if(*nrhs <= 1)
        {
            j = 1;
        L70: /* Solve U**T * x = b. */
            i__1 = j * b_dim1 + 1;
            z_div(&z__1, &b[j * b_dim1 + 1], &d__[1]);
            b[i__1].real = z__1.real;
            b[i__1].imag = z__1.imag; // , expr subst
            if(*n > 1)
            {
                i__1 = j * b_dim1 + 2;
                i__2 = j * b_dim1 + 2;
                i__3 = j * b_dim1 + 1;
                z__3.real = du[1].real * b[i__3].real - du[1].imag * b[i__3].imag;
                z__3.imag = du[1].real * b[i__3].imag + du[1].imag * b[i__3].real; // , expr subst
                z__2.real = b[i__2].real - z__3.real;
                z__2.imag = b[i__2].imag - z__3.imag; // , expr subst
                z_div(&z__1, &z__2, &d__[2]);
                b[i__1].real = z__1.real;
                b[i__1].imag = z__1.imag; // , expr subst
            }
            i__1 = *n;
            for(i__ = 3; i__ <= i__1; ++i__)
            {
                i__2 = i__ + j * b_dim1;
                i__3 = i__ + j * b_dim1;
                i__4 = i__ - 1;
                i__5 = i__ - 1 + j * b_dim1;
                z__4.real = du[i__4].real * b[i__5].real - du[i__4].imag * b[i__5].imag;
                z__4.imag = du[i__4].real * b[i__5].imag + du[i__4].imag * b[i__5].real; // , expr subst
                z__3.real = b[i__3].real - z__4.real;
                z__3.imag = b[i__3].imag - z__4.imag; // , expr subst
                i__6 = i__ - 2;
                i__7 = i__ - 2 + j * b_dim1;
                z__5.real = du2[i__6].real * b[i__7].real - du2[i__6].imag * b[i__7].imag;
                z__5.imag = du2[i__6].real * b[i__7].imag + du2[i__6].imag * b[i__7].real; // , expr subst
                z__2.real = z__3.real - z__5.real;
                z__2.imag = z__3.imag - z__5.imag; // , expr subst
                z_div(&z__1, &z__2, &d__[i__]);
                b[i__2].real = z__1.real;
                b[i__2].imag = z__1.imag; // , expr subst
                /* L80: */
            }
            /* Solve L**T * x = b. */
            for(i__ = *n - 1; i__ >= 1; --i__)
            {
                if(ipiv[i__] == i__)
                {
                    i__1 = i__ + j * b_dim1;
                    i__2 = i__ + j * b_dim1;
                    i__3 = i__;
                    i__4 = i__ + 1 + j * b_dim1;
                    z__2.real = dl[i__3].real * b[i__4].real - dl[i__3].imag * b[i__4].imag;
                    z__2.imag = dl[i__3].real * b[i__4].imag + dl[i__3].imag * b[i__4].real; // , expr subst
                    z__1.real = b[i__2].real - z__2.real;
                    z__1.imag = b[i__2].imag - z__2.imag; // , expr subst
                    b[i__1].real = z__1.real;
                    b[i__1].imag = z__1.imag; // , expr subst
                }
                else
                {
                    i__1 = i__ + 1 + j * b_dim1;
                    temp.real = b[i__1].real;
                    temp.imag = b[i__1].imag; // , expr subst
                    i__1 = i__ + 1 + j * b_dim1;
                    i__2 = i__ + j * b_dim1;
                    i__3 = i__;
                    z__2.real = dl[i__3].real * temp.real - dl[i__3].imag * temp.imag;
                    z__2.imag = dl[i__3].real * temp.imag + dl[i__3].imag * temp.real; // , expr subst
                    z__1.real = b[i__2].real - z__2.real;
                    z__1.imag = b[i__2].imag - z__2.imag; // , expr subst
                    b[i__1].real = z__1.real;
                    b[i__1].imag = z__1.imag; // , expr subst
                    i__1 = i__ + j * b_dim1;
                    b[i__1].real = temp.real;
                    b[i__1].imag = temp.imag; // , expr subst
                }
                /* L90: */
            }
            if(j < *nrhs)
            {
                ++j;
                goto L70;
            }
        }
        else
        {
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                /* Solve U**T * x = b. */
                i__2 = j * b_dim1 + 1;
                z_div(&z__1, &b[j * b_dim1 + 1], &d__[1]);
                b[i__2].real = z__1.real;
                b[i__2].imag = z__1.imag; // , expr subst
                if(*n > 1)
                {
                    i__2 = j * b_dim1 + 2;
                    i__3 = j * b_dim1 + 2;
                    i__4 = j * b_dim1 + 1;
                    z__3.real = du[1].real * b[i__4].real - du[1].imag * b[i__4].imag;
                    z__3.imag = du[1].real * b[i__4].imag + du[1].imag * b[i__4].real; // , expr subst
                    z__2.real = b[i__3].real - z__3.real;
                    z__2.imag = b[i__3].imag - z__3.imag; // , expr subst
                    z_div(&z__1, &z__2, &d__[2]);
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                }
                i__2 = *n;
                for(i__ = 3; i__ <= i__2; ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    i__4 = i__ + j * b_dim1;
                    i__5 = i__ - 1;
                    i__6 = i__ - 1 + j * b_dim1;
                    z__4.real = du[i__5].real * b[i__6].real - du[i__5].imag * b[i__6].imag;
                    z__4.imag = du[i__5].real * b[i__6].imag + du[i__5].imag * b[i__6].real; // , expr subst
                    z__3.real = b[i__4].real - z__4.real;
                    z__3.imag = b[i__4].imag - z__4.imag; // , expr subst
                    i__7 = i__ - 2;
                    i__8 = i__ - 2 + j * b_dim1;
                    z__5.real = du2[i__7].real * b[i__8].real - du2[i__7].imag * b[i__8].imag;
                    z__5.imag = du2[i__7].real * b[i__8].imag + du2[i__7].imag * b[i__8].real; // , expr subst
                    z__2.real = z__3.real - z__5.real;
                    z__2.imag = z__3.imag - z__5.imag; // , expr subst
                    z_div(&z__1, &z__2, &d__[i__]);
                    b[i__3].real = z__1.real;
                    b[i__3].imag = z__1.imag; // , expr subst
                    /* L100: */
                }
                /* Solve L**T * x = b. */
                for(i__ = *n - 1; i__ >= 1; --i__)
                {
                    if(ipiv[i__] == i__)
                    {
                        i__2 = i__ + j * b_dim1;
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__;
                        i__5 = i__ + 1 + j * b_dim1;
                        z__2.real = dl[i__4].real * b[i__5].real - dl[i__4].imag * b[i__5].imag;
                        z__2.imag = dl[i__4].real * b[i__5].imag + dl[i__4].imag * b[i__5].real; // , expr subst
                        z__1.real = b[i__3].real - z__2.real;
                        z__1.imag = b[i__3].imag - z__2.imag; // , expr subst
                        b[i__2].real = z__1.real;
                        b[i__2].imag = z__1.imag; // , expr subst
                    }
                    else
                    {
                        i__2 = i__ + 1 + j * b_dim1;
                        temp.real = b[i__2].real;
                        temp.imag = b[i__2].imag; // , expr subst
                        i__2 = i__ + 1 + j * b_dim1;
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__;
                        z__2.real = dl[i__4].real * temp.real - dl[i__4].imag * temp.imag;
                        z__2.imag = dl[i__4].real * temp.imag + dl[i__4].imag * temp.real; // , expr subst
                        z__1.real = b[i__3].real - z__2.real;
                        z__1.imag = b[i__3].imag - z__2.imag; // , expr subst
                        b[i__2].real = z__1.real;
                        b[i__2].imag = z__1.imag; // , expr subst
                        i__2 = i__ + j * b_dim1;
                        b[i__2].real = temp.real;
                        b[i__2].imag = temp.imag; // , expr subst
                    }
                    /* L110: */
                }
                /* L120: */
            }
        }
    }
    else
    {
        /* Solve A**H * X = B. */
        if(*nrhs <= 1)
        {
            j = 1;
        L130: /* Solve U**H * x = b. */
            i__1 = j * b_dim1 + 1;
            d_cnjg(&z__2, &d__[1]);
            z_div(&z__1, &b[j * b_dim1 + 1], &z__2);
            b[i__1].real = z__1.real;
            b[i__1].imag = z__1.imag; // , expr subst
            if(*n > 1)
            {
                i__1 = j * b_dim1 + 2;
                i__2 = j * b_dim1 + 2;
                d_cnjg(&z__4, &du[1]);
                i__3 = j * b_dim1 + 1;
                z__3.real = z__4.real * b[i__3].real - z__4.imag * b[i__3].imag;
                z__3.imag = z__4.real * b[i__3].imag + z__4.imag * b[i__3].real; // , expr subst
                z__2.real = b[i__2].real - z__3.real;
                z__2.imag = b[i__2].imag - z__3.imag; // , expr subst
                d_cnjg(&z__5, &d__[2]);
                z_div(&z__1, &z__2, &z__5);
                b[i__1].real = z__1.real;
                b[i__1].imag = z__1.imag; // , expr subst
            }
            i__1 = *n;
            for(i__ = 3; i__ <= i__1; ++i__)
            {
                i__2 = i__ + j * b_dim1;
                i__3 = i__ + j * b_dim1;
                d_cnjg(&z__5, &du[i__ - 1]);
                i__4 = i__ - 1 + j * b_dim1;
                z__4.real = z__5.real * b[i__4].real - z__5.imag * b[i__4].imag;
                z__4.imag = z__5.real * b[i__4].imag + z__5.imag * b[i__4].real; // , expr subst
                z__3.real = b[i__3].real - z__4.real;
                z__3.imag = b[i__3].imag - z__4.imag; // , expr subst
                d_cnjg(&z__7, &du2[i__ - 2]);
                i__5 = i__ - 2 + j * b_dim1;
                z__6.real = z__7.real * b[i__5].real - z__7.imag * b[i__5].imag;
                z__6.imag = z__7.real * b[i__5].imag + z__7.imag * b[i__5].real; // , expr subst
                z__2.real = z__3.real - z__6.real;
                z__2.imag = z__3.imag - z__6.imag; // , expr subst
                d_cnjg(&z__8, &d__[i__]);
                z_div(&z__1, &z__2, &z__8);
                b[i__2].real = z__1.real;
                b[i__2].imag = z__1.imag; // , expr subst
                /* L140: */
            }
            /* Solve L**H * x = b. */
            for(i__ = *n - 1; i__ >= 1; --i__)
            {
                if(ipiv[i__] == i__)
                {
                    i__1 = i__ + j * b_dim1;
                    i__2 = i__ + j * b_dim1;
                    d_cnjg(&z__3, &dl[i__]);
                    i__3 = i__ + 1 + j * b_dim1;
                    z__2.real = z__3.real * b[i__3].real - z__3.imag * b[i__3].imag;
                    z__2.imag = z__3.real * b[i__3].imag + z__3.imag * b[i__3].real; // , expr subst
                    z__1.real = b[i__2].real - z__2.real;
                    z__1.imag = b[i__2].imag - z__2.imag; // , expr subst
                    b[i__1].real = z__1.real;
                    b[i__1].imag = z__1.imag; // , expr subst
                }
                else
                {
                    i__1 = i__ + 1 + j * b_dim1;
                    temp.real = b[i__1].real;
                    temp.imag = b[i__1].imag; // , expr subst
                    i__1 = i__ + 1 + j * b_dim1;
                    i__2 = i__ + j * b_dim1;
                    d_cnjg(&z__3, &dl[i__]);
                    z__2.real = z__3.real * temp.real - z__3.imag * temp.imag;
                    z__2.imag = z__3.real * temp.imag + z__3.imag * temp.real; // , expr subst
                    z__1.real = b[i__2].real - z__2.real;
                    z__1.imag = b[i__2].imag - z__2.imag; // , expr subst
                    b[i__1].real = z__1.real;
                    b[i__1].imag = z__1.imag; // , expr subst
                    i__1 = i__ + j * b_dim1;
                    b[i__1].real = temp.real;
                    b[i__1].imag = temp.imag; // , expr subst
                }
                /* L150: */
            }
            if(j < *nrhs)
            {
                ++j;
                goto L130;
            }
        }
        else
        {
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                /* Solve U**H * x = b. */
                i__2 = j * b_dim1 + 1;
                d_cnjg(&z__2, &d__[1]);
                z_div(&z__1, &b[j * b_dim1 + 1], &z__2);
                b[i__2].real = z__1.real;
                b[i__2].imag = z__1.imag; // , expr subst
                if(*n > 1)
                {
                    i__2 = j * b_dim1 + 2;
                    i__3 = j * b_dim1 + 2;
                    d_cnjg(&z__4, &du[1]);
                    i__4 = j * b_dim1 + 1;
                    z__3.real = z__4.real * b[i__4].real - z__4.imag * b[i__4].imag;
                    z__3.imag = z__4.real * b[i__4].imag + z__4.imag * b[i__4].real; // , expr subst
                    z__2.real = b[i__3].real - z__3.real;
                    z__2.imag = b[i__3].imag - z__3.imag; // , expr subst
                    d_cnjg(&z__5, &d__[2]);
                    z_div(&z__1, &z__2, &z__5);
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                }
                i__2 = *n;
                for(i__ = 3; i__ <= i__2; ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    i__4 = i__ + j * b_dim1;
                    d_cnjg(&z__5, &du[i__ - 1]);
                    i__5 = i__ - 1 + j * b_dim1;
                    z__4.real = z__5.real * b[i__5].real - z__5.imag * b[i__5].imag;
                    z__4.imag = z__5.real * b[i__5].imag + z__5.imag * b[i__5].real; // , expr subst
                    z__3.real = b[i__4].real - z__4.real;
                    z__3.imag = b[i__4].imag - z__4.imag; // , expr subst
                    d_cnjg(&z__7, &du2[i__ - 2]);
                    i__6 = i__ - 2 + j * b_dim1;
                    z__6.real = z__7.real * b[i__6].real - z__7.imag * b[i__6].imag;
                    z__6.imag = z__7.real * b[i__6].imag + z__7.imag * b[i__6].real; // , expr subst
                    z__2.real = z__3.real - z__6.real;
                    z__2.imag = z__3.imag - z__6.imag; // , expr subst
                    d_cnjg(&z__8, &d__[i__]);
                    z_div(&z__1, &z__2, &z__8);
                    b[i__3].real = z__1.real;
                    b[i__3].imag = z__1.imag; // , expr subst
                    /* L160: */
                }
                /* Solve L**H * x = b. */
                for(i__ = *n - 1; i__ >= 1; --i__)
                {
                    if(ipiv[i__] == i__)
                    {
                        i__2 = i__ + j * b_dim1;
                        i__3 = i__ + j * b_dim1;
                        d_cnjg(&z__3, &dl[i__]);
                        i__4 = i__ + 1 + j * b_dim1;
                        z__2.real = z__3.real * b[i__4].real - z__3.imag * b[i__4].imag;
                        z__2.imag = z__3.real * b[i__4].imag + z__3.imag * b[i__4].real; // , expr subst
                        z__1.real = b[i__3].real - z__2.real;
                        z__1.imag = b[i__3].imag - z__2.imag; // , expr subst
                        b[i__2].real = z__1.real;
                        b[i__2].imag = z__1.imag; // , expr subst
                    }
                    else
                    {
                        i__2 = i__ + 1 + j * b_dim1;
                        temp.real = b[i__2].real;
                        temp.imag = b[i__2].imag; // , expr subst
                        i__2 = i__ + 1 + j * b_dim1;
                        i__3 = i__ + j * b_dim1;
                        d_cnjg(&z__3, &dl[i__]);
                        z__2.real = z__3.real * temp.real - z__3.imag * temp.imag;
                        z__2.imag = z__3.real * temp.imag + z__3.imag * temp.real; // , expr subst
                        z__1.real = b[i__3].real - z__2.real;
                        z__1.imag = b[i__3].imag - z__2.imag; // , expr subst
                        b[i__2].real = z__1.real;
                        b[i__2].imag = z__1.imag; // , expr subst
                        i__2 = i__ + j * b_dim1;
                        b[i__2].real = temp.real;
                        b[i__2].imag = temp.imag; // , expr subst
                    }
                    /* L170: */
                }
                /* L180: */
            }
        }
    }
    /* End of ZGTTS2 */
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
/* zgtts2_ */
