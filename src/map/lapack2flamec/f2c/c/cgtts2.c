/* ../netlib/cgtts2.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CGTTS2 solves a system of linear equations with a tridiagonal matrix using the LU factorization computed by sgttrf. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGTTS2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgtts2.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgtts2.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgtts2.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB ) */
/* .. Scalar Arguments .. */
/* INTEGER ITRANS, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGTTS2 solves one of the systems of equations */
/* > A * X = B, A**T * X = B, or A**H * X = B, */
/* > with a tridiagonal matrix A using the LU factorization computed */
/* > by CGTTRF. */
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
/* > DL is COMPLEX array, dimension (N-1) */
/* > The (n-1) multipliers that define the matrix L from the */
/* > LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is COMPLEX array, dimension (N) */
/* > The n diagonal elements of the upper triangular matrix U from */
/* > the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* > DU is COMPLEX array, dimension (N-1) */
/* > The (n-1) elements of the first super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* > DU2 is COMPLEX array, dimension (N-2) */
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
/* > B is COMPLEX array, dimension (LDB,NRHS) */
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
/* > \ingroup complexGTcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void cgtts2_(aocl_int_t *itrans, aocl_int_t *n, aocl_int_t *nrhs, scomplex *dl, scomplex *d__,
             scomplex *du, scomplex *du2, aocl_int_t *ipiv, scomplex *b, aocl_int_t *ldb)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgtts2(itrans, n, nrhs, dl, d__, du, du2, ipiv, b, ldb);
#else
    aocl_int64_t itrans_64 = *itrans;
    aocl_int64_t n_64 = *n;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t ldb_64 = *ldb;

    aocl_lapack_cgtts2(&itrans_64, &n_64, &nrhs_64, dl, d__, du, du2, ipiv, b, &ldb_64);
#endif
}

void aocl_lapack_cgtts2(aocl_int64_t *itrans, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *dl,
                        scomplex *d__, scomplex *du, scomplex *du2, aocl_int_t *ipiv, scomplex *b,
                        aocl_int64_t *ldb)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "cgtts2 inputs: itrans %lld, n %lld, nrhs %lld, ldb %lld", *itrans, *n,
             *nrhs, *ldb);
#else
    snprintf(buffer, 256, "cgtts2 inputs: itrans %d, n %d, nrhs %d, ldb %d", *itrans, *n, *nrhs,
             *ldb);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    scomplex q__1, q__2, q__3, q__4, q__5, q__6, q__7, q__8;
    /* Builtin functions */
    void c_div(scomplex *, scomplex *, scomplex *), r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    aocl_int64_t i__, j;
    scomplex temp;
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
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
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
                    q__2.real = dl[i__4].real * b[i__5].real - dl[i__4].imag * b[i__5].imag;
                    q__2.imag = dl[i__4].real * b[i__5].imag + dl[i__4].imag * b[i__5].real; // , expr subst
                    q__1.real = b[i__3].real - q__2.real;
                    q__1.imag = b[i__3].imag - q__2.imag; // , expr subst
                    b[i__2].real = q__1.real;
                    b[i__2].imag = q__1.imag; // , expr subst
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
                    q__2.real = dl[i__3].real * b[i__4].real - dl[i__3].imag * b[i__4].imag;
                    q__2.imag = dl[i__3].real * b[i__4].imag + dl[i__3].imag * b[i__4].real; // , expr subst
                    q__1.real = temp.real - q__2.real;
                    q__1.imag = temp.imag - q__2.imag; // , expr subst
                    b[i__2].real = q__1.real;
                    b[i__2].imag = q__1.imag; // , expr subst
                }
                /* L20: */
            }
            /* Solve U*x = b. */
            i__1 = *n + j * b_dim1;
            c_div(&q__1, &b[*n + j * b_dim1], &d__[*n]);
            b[i__1].real = q__1.real;
            b[i__1].imag = q__1.imag; // , expr subst
            if(*n > 1)
            {
                i__1 = *n - 1 + j * b_dim1;
                i__2 = *n - 1 + j * b_dim1;
                i__3 = *n - 1;
                i__4 = *n + j * b_dim1;
                q__3.real = du[i__3].real * b[i__4].real - du[i__3].imag * b[i__4].imag;
                q__3.imag = du[i__3].real * b[i__4].imag + du[i__3].imag * b[i__4].real; // , expr subst
                q__2.real = b[i__2].real - q__3.real;
                q__2.imag = b[i__2].imag - q__3.imag; // , expr subst
                c_div(&q__1, &q__2, &d__[*n - 1]);
                b[i__1].real = q__1.real;
                b[i__1].imag = q__1.imag; // , expr subst
            }
            for(i__ = *n - 2; i__ >= 1; --i__)
            {
                i__1 = i__ + j * b_dim1;
                i__2 = i__ + j * b_dim1;
                i__3 = i__;
                i__4 = i__ + 1 + j * b_dim1;
                q__4.real = du[i__3].real * b[i__4].real - du[i__3].imag * b[i__4].imag;
                q__4.imag = du[i__3].real * b[i__4].imag + du[i__3].imag * b[i__4].real; // , expr subst
                q__3.real = b[i__2].real - q__4.real;
                q__3.imag = b[i__2].imag - q__4.imag; // , expr subst
                i__5 = i__;
                i__6 = i__ + 2 + j * b_dim1;
                q__5.real = du2[i__5].real * b[i__6].real - du2[i__5].imag * b[i__6].imag;
                q__5.imag = du2[i__5].real * b[i__6].imag + du2[i__5].imag * b[i__6].real; // , expr subst
                q__2.real = q__3.real - q__5.real;
                q__2.imag = q__3.imag - q__5.imag; // , expr subst
                c_div(&q__1, &q__2, &d__[i__]);
                b[i__1].real = q__1.real;
                b[i__1].imag = q__1.imag; // , expr subst
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
                        q__2.real = dl[i__5].real * b[i__6].real - dl[i__5].imag * b[i__6].imag;
                        q__2.imag = dl[i__5].real * b[i__6].imag + dl[i__5].imag * b[i__6].real; // , expr subst
                        q__1.real = b[i__4].real - q__2.real;
                        q__1.imag = b[i__4].imag - q__2.imag; // , expr subst
                        b[i__3].real = q__1.real;
                        b[i__3].imag = q__1.imag; // , expr subst
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
                        q__2.real = dl[i__4].real * b[i__5].real - dl[i__4].imag * b[i__5].imag;
                        q__2.imag = dl[i__4].real * b[i__5].imag + dl[i__4].imag * b[i__5].real; // , expr subst
                        q__1.real = temp.real - q__2.real;
                        q__1.imag = temp.imag - q__2.imag; // , expr subst
                        b[i__3].real = q__1.real;
                        b[i__3].imag = q__1.imag; // , expr subst
                    }
                    /* L40: */
                }
                /* Solve U*x = b. */
                i__2 = *n + j * b_dim1;
                c_div(&q__1, &b[*n + j * b_dim1], &d__[*n]);
                b[i__2].real = q__1.real;
                b[i__2].imag = q__1.imag; // , expr subst
                if(*n > 1)
                {
                    i__2 = *n - 1 + j * b_dim1;
                    i__3 = *n - 1 + j * b_dim1;
                    i__4 = *n - 1;
                    i__5 = *n + j * b_dim1;
                    q__3.real = du[i__4].real * b[i__5].real - du[i__4].imag * b[i__5].imag;
                    q__3.imag = du[i__4].real * b[i__5].imag + du[i__4].imag * b[i__5].real; // , expr subst
                    q__2.real = b[i__3].real - q__3.real;
                    q__2.imag = b[i__3].imag - q__3.imag; // , expr subst
                    c_div(&q__1, &q__2, &d__[*n - 1]);
                    b[i__2].real = q__1.real;
                    b[i__2].imag = q__1.imag; // , expr subst
                }
                for(i__ = *n - 2; i__ >= 1; --i__)
                {
                    i__2 = i__ + j * b_dim1;
                    i__3 = i__ + j * b_dim1;
                    i__4 = i__;
                    i__5 = i__ + 1 + j * b_dim1;
                    q__4.real = du[i__4].real * b[i__5].real - du[i__4].imag * b[i__5].imag;
                    q__4.imag = du[i__4].real * b[i__5].imag + du[i__4].imag * b[i__5].real; // , expr subst
                    q__3.real = b[i__3].real - q__4.real;
                    q__3.imag = b[i__3].imag - q__4.imag; // , expr subst
                    i__6 = i__;
                    i__7 = i__ + 2 + j * b_dim1;
                    q__5.real = du2[i__6].real * b[i__7].real - du2[i__6].imag * b[i__7].imag;
                    q__5.imag = du2[i__6].real * b[i__7].imag + du2[i__6].imag * b[i__7].real; // , expr subst
                    q__2.real = q__3.real - q__5.real;
                    q__2.imag = q__3.imag - q__5.imag; // , expr subst
                    c_div(&q__1, &q__2, &d__[i__]);
                    b[i__2].real = q__1.real;
                    b[i__2].imag = q__1.imag; // , expr subst
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
            c_div(&q__1, &b[j * b_dim1 + 1], &d__[1]);
            b[i__1].real = q__1.real;
            b[i__1].imag = q__1.imag; // , expr subst
            if(*n > 1)
            {
                i__1 = j * b_dim1 + 2;
                i__2 = j * b_dim1 + 2;
                i__3 = j * b_dim1 + 1;
                q__3.real = du[1].real * b[i__3].real - du[1].imag * b[i__3].imag;
                q__3.imag = du[1].real * b[i__3].imag + du[1].imag * b[i__3].real; // , expr subst
                q__2.real = b[i__2].real - q__3.real;
                q__2.imag = b[i__2].imag - q__3.imag; // , expr subst
                c_div(&q__1, &q__2, &d__[2]);
                b[i__1].real = q__1.real;
                b[i__1].imag = q__1.imag; // , expr subst
            }
            i__1 = *n;
            for(i__ = 3; i__ <= i__1; ++i__)
            {
                i__2 = i__ + j * b_dim1;
                i__3 = i__ + j * b_dim1;
                i__4 = i__ - 1;
                i__5 = i__ - 1 + j * b_dim1;
                q__4.real = du[i__4].real * b[i__5].real - du[i__4].imag * b[i__5].imag;
                q__4.imag = du[i__4].real * b[i__5].imag + du[i__4].imag * b[i__5].real; // , expr subst
                q__3.real = b[i__3].real - q__4.real;
                q__3.imag = b[i__3].imag - q__4.imag; // , expr subst
                i__6 = i__ - 2;
                i__7 = i__ - 2 + j * b_dim1;
                q__5.real = du2[i__6].real * b[i__7].real - du2[i__6].imag * b[i__7].imag;
                q__5.imag = du2[i__6].real * b[i__7].imag + du2[i__6].imag * b[i__7].real; // , expr subst
                q__2.real = q__3.real - q__5.real;
                q__2.imag = q__3.imag - q__5.imag; // , expr subst
                c_div(&q__1, &q__2, &d__[i__]);
                b[i__2].real = q__1.real;
                b[i__2].imag = q__1.imag; // , expr subst
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
                    q__2.real = dl[i__3].real * b[i__4].real - dl[i__3].imag * b[i__4].imag;
                    q__2.imag = dl[i__3].real * b[i__4].imag + dl[i__3].imag * b[i__4].real; // , expr subst
                    q__1.real = b[i__2].real - q__2.real;
                    q__1.imag = b[i__2].imag - q__2.imag; // , expr subst
                    b[i__1].real = q__1.real;
                    b[i__1].imag = q__1.imag; // , expr subst
                }
                else
                {
                    i__1 = i__ + 1 + j * b_dim1;
                    temp.real = b[i__1].real;
                    temp.imag = b[i__1].imag; // , expr subst
                    i__1 = i__ + 1 + j * b_dim1;
                    i__2 = i__ + j * b_dim1;
                    i__3 = i__;
                    q__2.real = dl[i__3].real * temp.real - dl[i__3].imag * temp.imag;
                    q__2.imag = dl[i__3].real * temp.imag + dl[i__3].imag * temp.real; // , expr subst
                    q__1.real = b[i__2].real - q__2.real;
                    q__1.imag = b[i__2].imag - q__2.imag; // , expr subst
                    b[i__1].real = q__1.real;
                    b[i__1].imag = q__1.imag; // , expr subst
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
                c_div(&q__1, &b[j * b_dim1 + 1], &d__[1]);
                b[i__2].real = q__1.real;
                b[i__2].imag = q__1.imag; // , expr subst
                if(*n > 1)
                {
                    i__2 = j * b_dim1 + 2;
                    i__3 = j * b_dim1 + 2;
                    i__4 = j * b_dim1 + 1;
                    q__3.real = du[1].real * b[i__4].real - du[1].imag * b[i__4].imag;
                    q__3.imag = du[1].real * b[i__4].imag + du[1].imag * b[i__4].real; // , expr subst
                    q__2.real = b[i__3].real - q__3.real;
                    q__2.imag = b[i__3].imag - q__3.imag; // , expr subst
                    c_div(&q__1, &q__2, &d__[2]);
                    b[i__2].real = q__1.real;
                    b[i__2].imag = q__1.imag; // , expr subst
                }
                i__2 = *n;
                for(i__ = 3; i__ <= i__2; ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    i__4 = i__ + j * b_dim1;
                    i__5 = i__ - 1;
                    i__6 = i__ - 1 + j * b_dim1;
                    q__4.real = du[i__5].real * b[i__6].real - du[i__5].imag * b[i__6].imag;
                    q__4.imag = du[i__5].real * b[i__6].imag + du[i__5].imag * b[i__6].real; // , expr subst
                    q__3.real = b[i__4].real - q__4.real;
                    q__3.imag = b[i__4].imag - q__4.imag; // , expr subst
                    i__7 = i__ - 2;
                    i__8 = i__ - 2 + j * b_dim1;
                    q__5.real = du2[i__7].real * b[i__8].real - du2[i__7].imag * b[i__8].imag;
                    q__5.imag = du2[i__7].real * b[i__8].imag + du2[i__7].imag * b[i__8].real; // , expr subst
                    q__2.real = q__3.real - q__5.real;
                    q__2.imag = q__3.imag - q__5.imag; // , expr subst
                    c_div(&q__1, &q__2, &d__[i__]);
                    b[i__3].real = q__1.real;
                    b[i__3].imag = q__1.imag; // , expr subst
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
                        q__2.real = dl[i__4].real * b[i__5].real - dl[i__4].imag * b[i__5].imag;
                        q__2.imag = dl[i__4].real * b[i__5].imag + dl[i__4].imag * b[i__5].real; // , expr subst
                        q__1.real = b[i__3].real - q__2.real;
                        q__1.imag = b[i__3].imag - q__2.imag; // , expr subst
                        b[i__2].real = q__1.real;
                        b[i__2].imag = q__1.imag; // , expr subst
                    }
                    else
                    {
                        i__2 = i__ + 1 + j * b_dim1;
                        temp.real = b[i__2].real;
                        temp.imag = b[i__2].imag; // , expr subst
                        i__2 = i__ + 1 + j * b_dim1;
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__;
                        q__2.real = dl[i__4].real * temp.real - dl[i__4].imag * temp.imag;
                        q__2.imag = dl[i__4].real * temp.imag + dl[i__4].imag * temp.real; // , expr subst
                        q__1.real = b[i__3].real - q__2.real;
                        q__1.imag = b[i__3].imag - q__2.imag; // , expr subst
                        b[i__2].real = q__1.real;
                        b[i__2].imag = q__1.imag; // , expr subst
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
            r_cnjg(&q__2, &d__[1]);
            c_div(&q__1, &b[j * b_dim1 + 1], &q__2);
            b[i__1].real = q__1.real;
            b[i__1].imag = q__1.imag; // , expr subst
            if(*n > 1)
            {
                i__1 = j * b_dim1 + 2;
                i__2 = j * b_dim1 + 2;
                r_cnjg(&q__4, &du[1]);
                i__3 = j * b_dim1 + 1;
                q__3.real = q__4.real * b[i__3].real - q__4.imag * b[i__3].imag;
                q__3.imag = q__4.real * b[i__3].imag + q__4.imag * b[i__3].real; // , expr subst
                q__2.real = b[i__2].real - q__3.real;
                q__2.imag = b[i__2].imag - q__3.imag; // , expr subst
                r_cnjg(&q__5, &d__[2]);
                c_div(&q__1, &q__2, &q__5);
                b[i__1].real = q__1.real;
                b[i__1].imag = q__1.imag; // , expr subst
            }
            i__1 = *n;
            for(i__ = 3; i__ <= i__1; ++i__)
            {
                i__2 = i__ + j * b_dim1;
                i__3 = i__ + j * b_dim1;
                r_cnjg(&q__5, &du[i__ - 1]);
                i__4 = i__ - 1 + j * b_dim1;
                q__4.real = q__5.real * b[i__4].real - q__5.imag * b[i__4].imag;
                q__4.imag = q__5.real * b[i__4].imag + q__5.imag * b[i__4].real; // , expr subst
                q__3.real = b[i__3].real - q__4.real;
                q__3.imag = b[i__3].imag - q__4.imag; // , expr subst
                r_cnjg(&q__7, &du2[i__ - 2]);
                i__5 = i__ - 2 + j * b_dim1;
                q__6.real = q__7.real * b[i__5].real - q__7.imag * b[i__5].imag;
                q__6.imag = q__7.real * b[i__5].imag + q__7.imag * b[i__5].real; // , expr subst
                q__2.real = q__3.real - q__6.real;
                q__2.imag = q__3.imag - q__6.imag; // , expr subst
                r_cnjg(&q__8, &d__[i__]);
                c_div(&q__1, &q__2, &q__8);
                b[i__2].real = q__1.real;
                b[i__2].imag = q__1.imag; // , expr subst
                /* L140: */
            }
            /* Solve L**H * x = b. */
            for(i__ = *n - 1; i__ >= 1; --i__)
            {
                if(ipiv[i__] == i__)
                {
                    i__1 = i__ + j * b_dim1;
                    i__2 = i__ + j * b_dim1;
                    r_cnjg(&q__3, &dl[i__]);
                    i__3 = i__ + 1 + j * b_dim1;
                    q__2.real = q__3.real * b[i__3].real - q__3.imag * b[i__3].imag;
                    q__2.imag = q__3.real * b[i__3].imag + q__3.imag * b[i__3].real; // , expr subst
                    q__1.real = b[i__2].real - q__2.real;
                    q__1.imag = b[i__2].imag - q__2.imag; // , expr subst
                    b[i__1].real = q__1.real;
                    b[i__1].imag = q__1.imag; // , expr subst
                }
                else
                {
                    i__1 = i__ + 1 + j * b_dim1;
                    temp.real = b[i__1].real;
                    temp.imag = b[i__1].imag; // , expr subst
                    i__1 = i__ + 1 + j * b_dim1;
                    i__2 = i__ + j * b_dim1;
                    r_cnjg(&q__3, &dl[i__]);
                    q__2.real = q__3.real * temp.real - q__3.imag * temp.imag;
                    q__2.imag = q__3.real * temp.imag + q__3.imag * temp.real; // , expr subst
                    q__1.real = b[i__2].real - q__2.real;
                    q__1.imag = b[i__2].imag - q__2.imag; // , expr subst
                    b[i__1].real = q__1.real;
                    b[i__1].imag = q__1.imag; // , expr subst
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
                r_cnjg(&q__2, &d__[1]);
                c_div(&q__1, &b[j * b_dim1 + 1], &q__2);
                b[i__2].real = q__1.real;
                b[i__2].imag = q__1.imag; // , expr subst
                if(*n > 1)
                {
                    i__2 = j * b_dim1 + 2;
                    i__3 = j * b_dim1 + 2;
                    r_cnjg(&q__4, &du[1]);
                    i__4 = j * b_dim1 + 1;
                    q__3.real = q__4.real * b[i__4].real - q__4.imag * b[i__4].imag;
                    q__3.imag = q__4.real * b[i__4].imag + q__4.imag * b[i__4].real; // , expr subst
                    q__2.real = b[i__3].real - q__3.real;
                    q__2.imag = b[i__3].imag - q__3.imag; // , expr subst
                    r_cnjg(&q__5, &d__[2]);
                    c_div(&q__1, &q__2, &q__5);
                    b[i__2].real = q__1.real;
                    b[i__2].imag = q__1.imag; // , expr subst
                }
                i__2 = *n;
                for(i__ = 3; i__ <= i__2; ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    i__4 = i__ + j * b_dim1;
                    r_cnjg(&q__5, &du[i__ - 1]);
                    i__5 = i__ - 1 + j * b_dim1;
                    q__4.real = q__5.real * b[i__5].real - q__5.imag * b[i__5].imag;
                    q__4.imag = q__5.real * b[i__5].imag + q__5.imag * b[i__5].real; // , expr subst
                    q__3.real = b[i__4].real - q__4.real;
                    q__3.imag = b[i__4].imag - q__4.imag; // , expr subst
                    r_cnjg(&q__7, &du2[i__ - 2]);
                    i__6 = i__ - 2 + j * b_dim1;
                    q__6.real = q__7.real * b[i__6].real - q__7.imag * b[i__6].imag;
                    q__6.imag = q__7.real * b[i__6].imag + q__7.imag * b[i__6].real; // , expr subst
                    q__2.real = q__3.real - q__6.real;
                    q__2.imag = q__3.imag - q__6.imag; // , expr subst
                    r_cnjg(&q__8, &d__[i__]);
                    c_div(&q__1, &q__2, &q__8);
                    b[i__3].real = q__1.real;
                    b[i__3].imag = q__1.imag; // , expr subst
                    /* L160: */
                }
                /* Solve L**H * x = b. */
                for(i__ = *n - 1; i__ >= 1; --i__)
                {
                    if(ipiv[i__] == i__)
                    {
                        i__2 = i__ + j * b_dim1;
                        i__3 = i__ + j * b_dim1;
                        r_cnjg(&q__3, &dl[i__]);
                        i__4 = i__ + 1 + j * b_dim1;
                        q__2.real = q__3.real * b[i__4].real - q__3.imag * b[i__4].imag;
                        q__2.imag = q__3.real * b[i__4].imag + q__3.imag * b[i__4].real; // , expr subst
                        q__1.real = b[i__3].real - q__2.real;
                        q__1.imag = b[i__3].imag - q__2.imag; // , expr subst
                        b[i__2].real = q__1.real;
                        b[i__2].imag = q__1.imag; // , expr subst
                    }
                    else
                    {
                        i__2 = i__ + 1 + j * b_dim1;
                        temp.real = b[i__2].real;
                        temp.imag = b[i__2].imag; // , expr subst
                        i__2 = i__ + 1 + j * b_dim1;
                        i__3 = i__ + j * b_dim1;
                        r_cnjg(&q__3, &dl[i__]);
                        q__2.real = q__3.real * temp.real - q__3.imag * temp.imag;
                        q__2.imag = q__3.real * temp.imag + q__3.imag * temp.real; // , expr subst
                        q__1.real = b[i__3].real - q__2.real;
                        q__1.imag = b[i__3].imag - q__2.imag; // , expr subst
                        b[i__2].real = q__1.real;
                        b[i__2].imag = q__1.imag; // , expr subst
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
    /* End of CGTTS2 */
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
}
/* cgtts2_ */
