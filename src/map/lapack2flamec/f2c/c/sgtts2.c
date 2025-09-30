/* ../netlib/sgtts2.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SGTTS2 solves a system of linear equations with a tridiagonal matrix using the LU factorization computed by sgttrf. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGTTS2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgtts2.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgtts2.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgtts2.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB ) */
/* .. Scalar Arguments .. */
/* INTEGER ITRANS, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* REAL B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGTTS2 solves one of the systems of equations */
/* > A*X = B or A**T*X = B, */
/* > with a tridiagonal matrix A using the LU factorization computed */
/* > by SGTTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ITRANS */
/* > \verbatim */
/* > ITRANS is INTEGER */
/* > Specifies the form of the system of equations. */
/* > = 0: A * X = B (No transpose) */
/* > = 1: A**T* X = B (Transpose) */
/* > = 2: A**T* X = B (Conjugate transpose = Transpose) */
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
/* > DL is REAL array, dimension (N-1) */
/* > The (n-1) multipliers that define the matrix L from the */
/* > LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The n diagonal elements of the upper triangular matrix U from */
/* > the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* > DU is REAL array, dimension (N-1) */
/* > The (n-1) elements of the first super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* > DU2 is REAL array, dimension (N-2) */
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
/* > B is REAL array, dimension (LDB,NRHS) */
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
/* > \ingroup realGTcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void sgtts2_(aocl_int_t *itrans, aocl_int_t *n, aocl_int_t *nrhs, real *dl, real *d__, real *du,
             real *du2, aocl_int_t *ipiv, real *b, aocl_int_t *ldb)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgtts2(itrans, n, nrhs, dl, d__, du, du2, ipiv, b, ldb);
#else
    aocl_int64_t itrans_64 = *itrans;
    aocl_int64_t n_64 = *n;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t ldb_64 = *ldb;

    aocl_lapack_sgtts2(&itrans_64, &n_64, &nrhs_64, dl, d__, du, du2, ipiv, b, &ldb_64);
#endif
}

void aocl_lapack_sgtts2(aocl_int64_t *itrans, aocl_int64_t *n, aocl_int64_t *nrhs, real *dl,
                        real *d__, real *du, real *du2, aocl_int_t *ipiv, real *b,
                        aocl_int64_t *ldb)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgtts2 inputs: itrans %" FLA_IS ", n %" FLA_IS ", nrhs %" FLA_IS
                      ", ldb %" FLA_IS "",
                      *itrans, *n, *nrhs, *ldb);
    /* System generated locals */
    aocl_int64_t b_dim1, b_offset, i__1, i__2;
    /* Local variables */
    aocl_int64_t i__, j, ip;
    real temp;
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
                ip = ipiv[i__];
                temp = b[i__ + 1 - ip + i__ + j * b_dim1] - dl[i__] * b[ip + j * b_dim1];
                b[i__ + j * b_dim1] = b[ip + j * b_dim1];
                b[i__ + 1 + j * b_dim1] = temp;
                /* L20: */
            }
            /* Solve U*x = b. */
            b[*n + j * b_dim1] /= d__[*n];
            if(*n > 1)
            {
                b[*n - 1 + j * b_dim1]
                    = (b[*n - 1 + j * b_dim1] - du[*n - 1] * b[*n + j * b_dim1]) / d__[*n - 1];
            }
            for(i__ = *n - 2; i__ >= 1; --i__)
            {
                b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__] * b[i__ + 1 + j * b_dim1]
                                       - du2[i__] * b[i__ + 2 + j * b_dim1])
                                      / d__[i__];
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
                        b[i__ + 1 + j * b_dim1] -= dl[i__] * b[i__ + j * b_dim1];
                    }
                    else
                    {
                        temp = b[i__ + j * b_dim1];
                        b[i__ + j * b_dim1] = b[i__ + 1 + j * b_dim1];
                        b[i__ + 1 + j * b_dim1] = temp - dl[i__] * b[i__ + j * b_dim1];
                    }
                    /* L40: */
                }
                /* Solve U*x = b. */
                b[*n + j * b_dim1] /= d__[*n];
                if(*n > 1)
                {
                    b[*n - 1 + j * b_dim1]
                        = (b[*n - 1 + j * b_dim1] - du[*n - 1] * b[*n + j * b_dim1]) / d__[*n - 1];
                }
                for(i__ = *n - 2; i__ >= 1; --i__)
                {
                    b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__] * b[i__ + 1 + j * b_dim1]
                                           - du2[i__] * b[i__ + 2 + j * b_dim1])
                                          / d__[i__];
                    /* L50: */
                }
                /* L60: */
            }
        }
    }
    else
    {
        /* Solve A**T * X = B. */
        if(*nrhs <= 1)
        {
            /* Solve U**T*x = b. */
            j = 1;
        L70:
            b[j * b_dim1 + 1] /= d__[1];
            if(*n > 1)
            {
                b[j * b_dim1 + 2] = (b[j * b_dim1 + 2] - du[1] * b[j * b_dim1 + 1]) / d__[2];
            }
            i__1 = *n;
            for(i__ = 3; i__ <= i__1; ++i__)
            {
                b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__ - 1] * b[i__ - 1 + j * b_dim1]
                                       - du2[i__ - 2] * b[i__ - 2 + j * b_dim1])
                                      / d__[i__];
                /* L80: */
            }
            /* Solve L**T*x = b. */
            for(i__ = *n - 1; i__ >= 1; --i__)
            {
                ip = ipiv[i__];
                temp = b[i__ + j * b_dim1] - dl[i__] * b[i__ + 1 + j * b_dim1];
                b[i__ + j * b_dim1] = b[ip + j * b_dim1];
                b[ip + j * b_dim1] = temp;
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
                /* Solve U**T*x = b. */
                b[j * b_dim1 + 1] /= d__[1];
                if(*n > 1)
                {
                    b[j * b_dim1 + 2] = (b[j * b_dim1 + 2] - du[1] * b[j * b_dim1 + 1]) / d__[2];
                }
                i__2 = *n;
                for(i__ = 3; i__ <= i__2; ++i__)
                {
                    b[i__ + j * b_dim1]
                        = (b[i__ + j * b_dim1] - du[i__ - 1] * b[i__ - 1 + j * b_dim1]
                           - du2[i__ - 2] * b[i__ - 2 + j * b_dim1])
                          / d__[i__];
                    /* L100: */
                }
                for(i__ = *n - 1; i__ >= 1; --i__)
                {
                    if(ipiv[i__] == i__)
                    {
                        b[i__ + j * b_dim1] -= dl[i__] * b[i__ + 1 + j * b_dim1];
                    }
                    else
                    {
                        temp = b[i__ + 1 + j * b_dim1];
                        b[i__ + 1 + j * b_dim1] = b[i__ + j * b_dim1] - dl[i__] * temp;
                        b[i__ + j * b_dim1] = temp;
                    }
                    /* L110: */
                }
                /* L120: */
            }
        }
    }
    /* End of SGTTS2 */
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
/* sgtts2_ */
