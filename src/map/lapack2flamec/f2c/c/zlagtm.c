/* ../netlib/zlagtm.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZLAGTM performs a matrix-matrix product of the form C = αAB+βC, where A is a tridiagonal matr ix, B and C are rectangular matrices, and α and β are scalars, which may be 0, 1, or -1. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAGTM + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlagtm.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlagtm.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlagtm.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA, */
/* B, LDB ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER LDB, LDX, N, NRHS */
/* DOUBLE PRECISION ALPHA, BETA */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 B( LDB, * ), D( * ), DL( * ), DU( * ), */
/* $ X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAGTM performs a matrix-vector product of the form */
/* > */
/* > B := alpha * A * X + beta * B */
/* > */
/* > where A is a tridiagonal matrix of order N, B and X are N by NRHS */
/* > matrices, and alpha and beta are real scalars, each of which may be */
/* > 0., 1., or -1. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the operation applied to A. */
/* > = 'N': No transpose, B := alpha * A * X + beta * B */
/* > = 'T': Transpose, B := alpha * A**T * X + beta * B */
/* > = 'C': Conjugate transpose, B := alpha * A**H * X + beta * B */
/* > \endverbatim */
/* > */
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
/* > of the matrices X and B. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* > ALPHA is DOUBLE PRECISION */
/* > The scalar alpha. ALPHA must be 0., 1., or -1.;
otherwise, */
/* > it is assumed to be 0. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* > DL is COMPLEX*16 array, dimension (N-1) */
/* > The (n-1) sub-diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is COMPLEX*16 array, dimension (N) */
/* > The diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* > DU is COMPLEX*16 array, dimension (N-1) */
/* > The (n-1) super-diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* > X is COMPLEX*16 array, dimension (LDX,NRHS) */
/* > The N by NRHS matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX is INTEGER */
/* > The leading dimension of the array X. LDX >= fla_max(N,1). */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* > BETA is DOUBLE PRECISION */
/* > The scalar beta. BETA must be 0., 1., or -1.;
otherwise, */
/* > it is assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* > On entry, the N by NRHS matrix B. */
/* > On exit, B is overwritten by the matrix expression */
/* > B := alpha * A * X + beta * B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(N,1). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zlagtm_(char *trans, aocl_int_t *n, aocl_int_t *nrhs, doublereal *alpha, dcomplex *dl,
             dcomplex *d__, dcomplex *du, dcomplex *x, aocl_int_t *ldx,
             doublereal *beta, dcomplex *b, aocl_int_t *ldb)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zlagtm(trans, n, nrhs, alpha, dl, d__, du, x, ldx, beta, b, ldb);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t ldx_64 = *ldx;
    aocl_int64_t ldb_64 = *ldb;

    aocl_lapack_zlagtm(trans, &n_64, &nrhs_64, alpha, dl, d__, du, x, &ldx_64, beta, b, &ldb_64);
#endif
}

void aocl_lapack_zlagtm(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, doublereal *alpha,
                        dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *x,
                        aocl_int64_t *ldx, doublereal *beta, dcomplex *b, aocl_int64_t *ldb)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlagtm inputs: trans %c, n %" FLA_IS ", nrhs %" FLA_IS ", ldx %" FLA_IS
                      ", ldb %" FLA_IS "",
                      *trans, *n, *nrhs, *ldx, *ldb);
    /* System generated locals */
    aocl_int64_t b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8,
        i__9, i__10;
    dcomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8, z__9;
    /* Builtin functions */
    void d_cnjg(dcomplex *, dcomplex *);
    /* Local variables */
    aocl_int64_t i__, j;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
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
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --dl;
    --d__;
    --du;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Multiply B by BETA if BETA.NE.1. */
    if(*beta == 0.)
    {
        i__1 = *nrhs;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *n;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                i__3 = i__ + j * b_dim1;
                b[i__3].real = 0.;
                b[i__3].imag = 0.; // , expr subst
                /* L10: */
            }
            /* L20: */
        }
    }
    else if(*beta == -1.)
    {
        i__1 = *nrhs;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *n;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                i__3 = i__ + j * b_dim1;
                i__4 = i__ + j * b_dim1;
                z__1.real = -b[i__4].real;
                z__1.imag = -b[i__4].imag; // , expr subst
                b[i__3].real = z__1.real;
                b[i__3].imag = z__1.imag; // , expr subst
                /* L30: */
            }
            /* L40: */
        }
    }
    if(*alpha == 1.)
    {
        if(lsame_(trans, "N", 1, 1))
        {
            /* Compute B := B + A*X */
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                if(*n == 1)
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    z__2.real = d__[1].real * x[i__4].real - d__[1].imag * x[i__4].imag;
                    z__2.imag = d__[1].real * x[i__4].imag + d__[1].imag * x[i__4].real; // , expr subst
                    z__1.real = b[i__3].real + z__2.real;
                    z__1.imag = b[i__3].imag + z__2.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                }
                else
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    z__3.real = d__[1].real * x[i__4].real - d__[1].imag * x[i__4].imag;
                    z__3.imag = d__[1].real * x[i__4].imag + d__[1].imag * x[i__4].real; // , expr subst
                    z__2.real = b[i__3].real + z__3.real;
                    z__2.imag = b[i__3].imag + z__3.imag; // , expr subst
                    i__5 = j * x_dim1 + 2;
                    z__4.real = du[1].real * x[i__5].real - du[1].imag * x[i__5].imag;
                    z__4.imag = du[1].real * x[i__5].imag + du[1].imag * x[i__5].real; // , expr subst
                    z__1.real = z__2.real + z__4.real;
                    z__1.imag = z__2.imag + z__4.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                    i__2 = *n + j * b_dim1;
                    i__3 = *n + j * b_dim1;
                    i__4 = *n - 1;
                    i__5 = *n - 1 + j * x_dim1;
                    z__3.real = dl[i__4].real * x[i__5].real - dl[i__4].imag * x[i__5].imag;
                    z__3.imag = dl[i__4].real * x[i__5].imag + dl[i__4].imag * x[i__5].real; // , expr subst
                    z__2.real = b[i__3].real + z__3.real;
                    z__2.imag = b[i__3].imag + z__3.imag; // , expr subst
                    i__6 = *n;
                    i__7 = *n + j * x_dim1;
                    z__4.real = d__[i__6].real * x[i__7].real - d__[i__6].imag * x[i__7].imag;
                    z__4.imag = d__[i__6].real * x[i__7].imag + d__[i__6].imag * x[i__7].real; // , expr subst
                    z__1.real = z__2.real + z__4.real;
                    z__1.imag = z__2.imag + z__4.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                    i__2 = *n - 1;
                    for(i__ = 2; i__ <= i__2; ++i__)
                    {
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        i__5 = i__ - 1;
                        i__6 = i__ - 1 + j * x_dim1;
                        z__4.real = dl[i__5].real * x[i__6].real - dl[i__5].imag * x[i__6].imag;
                        z__4.imag = dl[i__5].real * x[i__6].imag + dl[i__5].imag * x[i__6].real; // , expr subst
                        z__3.real = b[i__4].real + z__4.real;
                        z__3.imag = b[i__4].imag + z__4.imag; // , expr subst
                        i__7 = i__;
                        i__8 = i__ + j * x_dim1;
                        z__5.real = d__[i__7].real * x[i__8].real - d__[i__7].imag * x[i__8].imag;
                        z__5.imag = d__[i__7].real * x[i__8].imag + d__[i__7].imag * x[i__8].real; // , expr subst
                        z__2.real = z__3.real + z__5.real;
                        z__2.imag = z__3.imag + z__5.imag; // , expr subst
                        i__9 = i__;
                        i__10 = i__ + 1 + j * x_dim1;
                        z__6.real = du[i__9].real * x[i__10].real - du[i__9].imag * x[i__10].imag;
                        z__6.imag = du[i__9].real * x[i__10].imag + du[i__9].imag * x[i__10].real; // , expr subst
                        z__1.real = z__2.real + z__6.real;
                        z__1.imag = z__2.imag + z__6.imag; // , expr subst
                        b[i__3].real = z__1.real;
                        b[i__3].imag = z__1.imag; // , expr subst
                        /* L50: */
                    }
                }
                /* L60: */
            }
        }
        else if(lsame_(trans, "T", 1, 1))
        {
            /* Compute B := B + A**T * X */
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                if(*n == 1)
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    z__2.real = d__[1].real * x[i__4].real - d__[1].imag * x[i__4].imag;
                    z__2.imag = d__[1].real * x[i__4].imag + d__[1].imag * x[i__4].real; // , expr subst
                    z__1.real = b[i__3].real + z__2.real;
                    z__1.imag = b[i__3].imag + z__2.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                }
                else
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    z__3.real = d__[1].real * x[i__4].real - d__[1].imag * x[i__4].imag;
                    z__3.imag = d__[1].real * x[i__4].imag + d__[1].imag * x[i__4].real; // , expr subst
                    z__2.real = b[i__3].real + z__3.real;
                    z__2.imag = b[i__3].imag + z__3.imag; // , expr subst
                    i__5 = j * x_dim1 + 2;
                    z__4.real = dl[1].real * x[i__5].real - dl[1].imag * x[i__5].imag;
                    z__4.imag = dl[1].real * x[i__5].imag + dl[1].imag * x[i__5].real; // , expr subst
                    z__1.real = z__2.real + z__4.real;
                    z__1.imag = z__2.imag + z__4.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                    i__2 = *n + j * b_dim1;
                    i__3 = *n + j * b_dim1;
                    i__4 = *n - 1;
                    i__5 = *n - 1 + j * x_dim1;
                    z__3.real = du[i__4].real * x[i__5].real - du[i__4].imag * x[i__5].imag;
                    z__3.imag = du[i__4].real * x[i__5].imag + du[i__4].imag * x[i__5].real; // , expr subst
                    z__2.real = b[i__3].real + z__3.real;
                    z__2.imag = b[i__3].imag + z__3.imag; // , expr subst
                    i__6 = *n;
                    i__7 = *n + j * x_dim1;
                    z__4.real = d__[i__6].real * x[i__7].real - d__[i__6].imag * x[i__7].imag;
                    z__4.imag = d__[i__6].real * x[i__7].imag + d__[i__6].imag * x[i__7].real; // , expr subst
                    z__1.real = z__2.real + z__4.real;
                    z__1.imag = z__2.imag + z__4.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                    i__2 = *n - 1;
                    for(i__ = 2; i__ <= i__2; ++i__)
                    {
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        i__5 = i__ - 1;
                        i__6 = i__ - 1 + j * x_dim1;
                        z__4.real = du[i__5].real * x[i__6].real - du[i__5].imag * x[i__6].imag;
                        z__4.imag = du[i__5].real * x[i__6].imag + du[i__5].imag * x[i__6].real; // , expr subst
                        z__3.real = b[i__4].real + z__4.real;
                        z__3.imag = b[i__4].imag + z__4.imag; // , expr subst
                        i__7 = i__;
                        i__8 = i__ + j * x_dim1;
                        z__5.real = d__[i__7].real * x[i__8].real - d__[i__7].imag * x[i__8].imag;
                        z__5.imag = d__[i__7].real * x[i__8].imag + d__[i__7].imag * x[i__8].real; // , expr subst
                        z__2.real = z__3.real + z__5.real;
                        z__2.imag = z__3.imag + z__5.imag; // , expr subst
                        i__9 = i__;
                        i__10 = i__ + 1 + j * x_dim1;
                        z__6.real = dl[i__9].real * x[i__10].real - dl[i__9].imag * x[i__10].imag;
                        z__6.imag = dl[i__9].real * x[i__10].imag + dl[i__9].imag * x[i__10].real; // , expr subst
                        z__1.real = z__2.real + z__6.real;
                        z__1.imag = z__2.imag + z__6.imag; // , expr subst
                        b[i__3].real = z__1.real;
                        b[i__3].imag = z__1.imag; // , expr subst
                        /* L70: */
                    }
                }
                /* L80: */
            }
        }
        else if(lsame_(trans, "C", 1, 1))
        {
            /* Compute B := B + A**H * X */
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                if(*n == 1)
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    d_cnjg(&z__3, &d__[1]);
                    i__4 = j * x_dim1 + 1;
                    z__2.real = z__3.real * x[i__4].real - z__3.imag * x[i__4].imag;
                    z__2.imag = z__3.real * x[i__4].imag + z__3.imag * x[i__4].real; // , expr subst
                    z__1.real = b[i__3].real + z__2.real;
                    z__1.imag = b[i__3].imag + z__2.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                }
                else
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    d_cnjg(&z__4, &d__[1]);
                    i__4 = j * x_dim1 + 1;
                    z__3.real = z__4.real * x[i__4].real - z__4.imag * x[i__4].imag;
                    z__3.imag = z__4.real * x[i__4].imag + z__4.imag * x[i__4].real; // , expr subst
                    z__2.real = b[i__3].real + z__3.real;
                    z__2.imag = b[i__3].imag + z__3.imag; // , expr subst
                    d_cnjg(&z__6, &dl[1]);
                    i__5 = j * x_dim1 + 2;
                    z__5.real = z__6.real * x[i__5].real - z__6.imag * x[i__5].imag;
                    z__5.imag = z__6.real * x[i__5].imag + z__6.imag * x[i__5].real; // , expr subst
                    z__1.real = z__2.real + z__5.real;
                    z__1.imag = z__2.imag + z__5.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                    i__2 = *n + j * b_dim1;
                    i__3 = *n + j * b_dim1;
                    d_cnjg(&z__4, &du[*n - 1]);
                    i__4 = *n - 1 + j * x_dim1;
                    z__3.real = z__4.real * x[i__4].real - z__4.imag * x[i__4].imag;
                    z__3.imag = z__4.real * x[i__4].imag + z__4.imag * x[i__4].real; // , expr subst
                    z__2.real = b[i__3].real + z__3.real;
                    z__2.imag = b[i__3].imag + z__3.imag; // , expr subst
                    d_cnjg(&z__6, &d__[*n]);
                    i__5 = *n + j * x_dim1;
                    z__5.real = z__6.real * x[i__5].real - z__6.imag * x[i__5].imag;
                    z__5.imag = z__6.real * x[i__5].imag + z__6.imag * x[i__5].real; // , expr subst
                    z__1.real = z__2.real + z__5.real;
                    z__1.imag = z__2.imag + z__5.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                    i__2 = *n - 1;
                    for(i__ = 2; i__ <= i__2; ++i__)
                    {
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        d_cnjg(&z__5, &du[i__ - 1]);
                        i__5 = i__ - 1 + j * x_dim1;
                        z__4.real = z__5.real * x[i__5].real - z__5.imag * x[i__5].imag;
                        z__4.imag = z__5.real * x[i__5].imag + z__5.imag * x[i__5].real; // , expr subst
                        z__3.real = b[i__4].real + z__4.real;
                        z__3.imag = b[i__4].imag + z__4.imag; // , expr subst
                        d_cnjg(&z__7, &d__[i__]);
                        i__6 = i__ + j * x_dim1;
                        z__6.real = z__7.real * x[i__6].real - z__7.imag * x[i__6].imag;
                        z__6.imag = z__7.real * x[i__6].imag + z__7.imag * x[i__6].real; // , expr subst
                        z__2.real = z__3.real + z__6.real;
                        z__2.imag = z__3.imag + z__6.imag; // , expr subst
                        d_cnjg(&z__9, &dl[i__]);
                        i__7 = i__ + 1 + j * x_dim1;
                        z__8.real = z__9.real * x[i__7].real - z__9.imag * x[i__7].imag;
                        z__8.imag = z__9.real * x[i__7].imag + z__9.imag * x[i__7].real; // , expr subst
                        z__1.real = z__2.real + z__8.real;
                        z__1.imag = z__2.imag + z__8.imag; // , expr subst
                        b[i__3].real = z__1.real;
                        b[i__3].imag = z__1.imag; // , expr subst
                        /* L90: */
                    }
                }
                /* L100: */
            }
        }
    }
    else if(*alpha == -1.)
    {
        if(lsame_(trans, "N", 1, 1))
        {
            /* Compute B := B - A*X */
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                if(*n == 1)
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    z__2.real = d__[1].real * x[i__4].real - d__[1].imag * x[i__4].imag;
                    z__2.imag = d__[1].real * x[i__4].imag + d__[1].imag * x[i__4].real; // , expr subst
                    z__1.real = b[i__3].real - z__2.real;
                    z__1.imag = b[i__3].imag - z__2.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                }
                else
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    z__3.real = d__[1].real * x[i__4].real - d__[1].imag * x[i__4].imag;
                    z__3.imag = d__[1].real * x[i__4].imag + d__[1].imag * x[i__4].real; // , expr subst
                    z__2.real = b[i__3].real - z__3.real;
                    z__2.imag = b[i__3].imag - z__3.imag; // , expr subst
                    i__5 = j * x_dim1 + 2;
                    z__4.real = du[1].real * x[i__5].real - du[1].imag * x[i__5].imag;
                    z__4.imag = du[1].real * x[i__5].imag + du[1].imag * x[i__5].real; // , expr subst
                    z__1.real = z__2.real - z__4.real;
                    z__1.imag = z__2.imag - z__4.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                    i__2 = *n + j * b_dim1;
                    i__3 = *n + j * b_dim1;
                    i__4 = *n - 1;
                    i__5 = *n - 1 + j * x_dim1;
                    z__3.real = dl[i__4].real * x[i__5].real - dl[i__4].imag * x[i__5].imag;
                    z__3.imag = dl[i__4].real * x[i__5].imag + dl[i__4].imag * x[i__5].real; // , expr subst
                    z__2.real = b[i__3].real - z__3.real;
                    z__2.imag = b[i__3].imag - z__3.imag; // , expr subst
                    i__6 = *n;
                    i__7 = *n + j * x_dim1;
                    z__4.real = d__[i__6].real * x[i__7].real - d__[i__6].imag * x[i__7].imag;
                    z__4.imag = d__[i__6].real * x[i__7].imag + d__[i__6].imag * x[i__7].real; // , expr subst
                    z__1.real = z__2.real - z__4.real;
                    z__1.imag = z__2.imag - z__4.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                    i__2 = *n - 1;
                    for(i__ = 2; i__ <= i__2; ++i__)
                    {
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        i__5 = i__ - 1;
                        i__6 = i__ - 1 + j * x_dim1;
                        z__4.real = dl[i__5].real * x[i__6].real - dl[i__5].imag * x[i__6].imag;
                        z__4.imag = dl[i__5].real * x[i__6].imag + dl[i__5].imag * x[i__6].real; // , expr subst
                        z__3.real = b[i__4].real - z__4.real;
                        z__3.imag = b[i__4].imag - z__4.imag; // , expr subst
                        i__7 = i__;
                        i__8 = i__ + j * x_dim1;
                        z__5.real = d__[i__7].real * x[i__8].real - d__[i__7].imag * x[i__8].imag;
                        z__5.imag = d__[i__7].real * x[i__8].imag + d__[i__7].imag * x[i__8].real; // , expr subst
                        z__2.real = z__3.real - z__5.real;
                        z__2.imag = z__3.imag - z__5.imag; // , expr subst
                        i__9 = i__;
                        i__10 = i__ + 1 + j * x_dim1;
                        z__6.real = du[i__9].real * x[i__10].real - du[i__9].imag * x[i__10].imag;
                        z__6.imag = du[i__9].real * x[i__10].imag + du[i__9].imag * x[i__10].real; // , expr subst
                        z__1.real = z__2.real - z__6.real;
                        z__1.imag = z__2.imag - z__6.imag; // , expr subst
                        b[i__3].real = z__1.real;
                        b[i__3].imag = z__1.imag; // , expr subst
                        /* L110: */
                    }
                }
                /* L120: */
            }
        }
        else if(lsame_(trans, "T", 1, 1))
        {
            /* Compute B := B - A**T *X */
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                if(*n == 1)
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    z__2.real = d__[1].real * x[i__4].real - d__[1].imag * x[i__4].imag;
                    z__2.imag = d__[1].real * x[i__4].imag + d__[1].imag * x[i__4].real; // , expr subst
                    z__1.real = b[i__3].real - z__2.real;
                    z__1.imag = b[i__3].imag - z__2.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                }
                else
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    z__3.real = d__[1].real * x[i__4].real - d__[1].imag * x[i__4].imag;
                    z__3.imag = d__[1].real * x[i__4].imag + d__[1].imag * x[i__4].real; // , expr subst
                    z__2.real = b[i__3].real - z__3.real;
                    z__2.imag = b[i__3].imag - z__3.imag; // , expr subst
                    i__5 = j * x_dim1 + 2;
                    z__4.real = dl[1].real * x[i__5].real - dl[1].imag * x[i__5].imag;
                    z__4.imag = dl[1].real * x[i__5].imag + dl[1].imag * x[i__5].real; // , expr subst
                    z__1.real = z__2.real - z__4.real;
                    z__1.imag = z__2.imag - z__4.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                    i__2 = *n + j * b_dim1;
                    i__3 = *n + j * b_dim1;
                    i__4 = *n - 1;
                    i__5 = *n - 1 + j * x_dim1;
                    z__3.real = du[i__4].real * x[i__5].real - du[i__4].imag * x[i__5].imag;
                    z__3.imag = du[i__4].real * x[i__5].imag + du[i__4].imag * x[i__5].real; // , expr subst
                    z__2.real = b[i__3].real - z__3.real;
                    z__2.imag = b[i__3].imag - z__3.imag; // , expr subst
                    i__6 = *n;
                    i__7 = *n + j * x_dim1;
                    z__4.real = d__[i__6].real * x[i__7].real - d__[i__6].imag * x[i__7].imag;
                    z__4.imag = d__[i__6].real * x[i__7].imag + d__[i__6].imag * x[i__7].real; // , expr subst
                    z__1.real = z__2.real - z__4.real;
                    z__1.imag = z__2.imag - z__4.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                    i__2 = *n - 1;
                    for(i__ = 2; i__ <= i__2; ++i__)
                    {
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        i__5 = i__ - 1;
                        i__6 = i__ - 1 + j * x_dim1;
                        z__4.real = du[i__5].real * x[i__6].real - du[i__5].imag * x[i__6].imag;
                        z__4.imag = du[i__5].real * x[i__6].imag + du[i__5].imag * x[i__6].real; // , expr subst
                        z__3.real = b[i__4].real - z__4.real;
                        z__3.imag = b[i__4].imag - z__4.imag; // , expr subst
                        i__7 = i__;
                        i__8 = i__ + j * x_dim1;
                        z__5.real = d__[i__7].real * x[i__8].real - d__[i__7].imag * x[i__8].imag;
                        z__5.imag = d__[i__7].real * x[i__8].imag + d__[i__7].imag * x[i__8].real; // , expr subst
                        z__2.real = z__3.real - z__5.real;
                        z__2.imag = z__3.imag - z__5.imag; // , expr subst
                        i__9 = i__;
                        i__10 = i__ + 1 + j * x_dim1;
                        z__6.real = dl[i__9].real * x[i__10].real - dl[i__9].imag * x[i__10].imag;
                        z__6.imag = dl[i__9].real * x[i__10].imag + dl[i__9].imag * x[i__10].real; // , expr subst
                        z__1.real = z__2.real - z__6.real;
                        z__1.imag = z__2.imag - z__6.imag; // , expr subst
                        b[i__3].real = z__1.real;
                        b[i__3].imag = z__1.imag; // , expr subst
                        /* L130: */
                    }
                }
                /* L140: */
            }
        }
        else if(lsame_(trans, "C", 1, 1))
        {
            /* Compute B := B - A**H *X */
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                if(*n == 1)
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    d_cnjg(&z__3, &d__[1]);
                    i__4 = j * x_dim1 + 1;
                    z__2.real = z__3.real * x[i__4].real - z__3.imag * x[i__4].imag;
                    z__2.imag = z__3.real * x[i__4].imag + z__3.imag * x[i__4].real; // , expr subst
                    z__1.real = b[i__3].real - z__2.real;
                    z__1.imag = b[i__3].imag - z__2.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                }
                else
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    d_cnjg(&z__4, &d__[1]);
                    i__4 = j * x_dim1 + 1;
                    z__3.real = z__4.real * x[i__4].real - z__4.imag * x[i__4].imag;
                    z__3.imag = z__4.real * x[i__4].imag + z__4.imag * x[i__4].real; // , expr subst
                    z__2.real = b[i__3].real - z__3.real;
                    z__2.imag = b[i__3].imag - z__3.imag; // , expr subst
                    d_cnjg(&z__6, &dl[1]);
                    i__5 = j * x_dim1 + 2;
                    z__5.real = z__6.real * x[i__5].real - z__6.imag * x[i__5].imag;
                    z__5.imag = z__6.real * x[i__5].imag + z__6.imag * x[i__5].real; // , expr subst
                    z__1.real = z__2.real - z__5.real;
                    z__1.imag = z__2.imag - z__5.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                    i__2 = *n + j * b_dim1;
                    i__3 = *n + j * b_dim1;
                    d_cnjg(&z__4, &du[*n - 1]);
                    i__4 = *n - 1 + j * x_dim1;
                    z__3.real = z__4.real * x[i__4].real - z__4.imag * x[i__4].imag;
                    z__3.imag = z__4.real * x[i__4].imag + z__4.imag * x[i__4].real; // , expr subst
                    z__2.real = b[i__3].real - z__3.real;
                    z__2.imag = b[i__3].imag - z__3.imag; // , expr subst
                    d_cnjg(&z__6, &d__[*n]);
                    i__5 = *n + j * x_dim1;
                    z__5.real = z__6.real * x[i__5].real - z__6.imag * x[i__5].imag;
                    z__5.imag = z__6.real * x[i__5].imag + z__6.imag * x[i__5].real; // , expr subst
                    z__1.real = z__2.real - z__5.real;
                    z__1.imag = z__2.imag - z__5.imag; // , expr subst
                    b[i__2].real = z__1.real;
                    b[i__2].imag = z__1.imag; // , expr subst
                    i__2 = *n - 1;
                    for(i__ = 2; i__ <= i__2; ++i__)
                    {
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        d_cnjg(&z__5, &du[i__ - 1]);
                        i__5 = i__ - 1 + j * x_dim1;
                        z__4.real = z__5.real * x[i__5].real - z__5.imag * x[i__5].imag;
                        z__4.imag = z__5.real * x[i__5].imag + z__5.imag * x[i__5].real; // , expr subst
                        z__3.real = b[i__4].real - z__4.real;
                        z__3.imag = b[i__4].imag - z__4.imag; // , expr subst
                        d_cnjg(&z__7, &d__[i__]);
                        i__6 = i__ + j * x_dim1;
                        z__6.real = z__7.real * x[i__6].real - z__7.imag * x[i__6].imag;
                        z__6.imag = z__7.real * x[i__6].imag + z__7.imag * x[i__6].real; // , expr subst
                        z__2.real = z__3.real - z__6.real;
                        z__2.imag = z__3.imag - z__6.imag; // , expr subst
                        d_cnjg(&z__9, &dl[i__]);
                        i__7 = i__ + 1 + j * x_dim1;
                        z__8.real = z__9.real * x[i__7].real - z__9.imag * x[i__7].imag;
                        z__8.imag = z__9.real * x[i__7].imag + z__9.imag * x[i__7].real; // , expr subst
                        z__1.real = z__2.real - z__8.real;
                        z__1.imag = z__2.imag - z__8.imag; // , expr subst
                        b[i__3].real = z__1.real;
                        b[i__3].imag = z__1.imag; // , expr subst
                        /* L150: */
                    }
                }
                /* L160: */
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZLAGTM */
}
/* zlagtm_ */
