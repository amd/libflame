/* ../netlib/sptts2.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SPTTS2 solves a tridiagonal system of the form AX=B using the L D LH factorization computed by spttrf. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SPTTS2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sptts2.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sptts2.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sptts2.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SPTTS2( N, NRHS, D, E, B, LDB ) */
/* .. Scalar Arguments .. */
/* INTEGER LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* REAL B( LDB, * ), D( * ), E( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPTTS2 solves a tridiagonal system of the form */
/* > A * X = B */
/* > using the L*D*L**T factorization of A computed by SPTTRF. D is a */
/* > diagonal matrix specified in the vector D, L is a unit bidiagonal */
/* > matrix whose subdiagonal is specified in the vector E, and X and B */
/* > are N by NRHS matrices. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the tridiagonal matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The n diagonal elements of the diagonal matrix D from the */
/* > L*D*L**T factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is REAL array, dimension (N-1) */
/* > The (n-1) subdiagonal elements of the unit bidiagonal factor */
/* > L from the L*D*L**T factorization of A. E can also be regarded */
/* > as the superdiagonal of the unit bidiagonal factor U from the */
/* > factorization A = U**T*D*U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB,NRHS) */
/* > On entry, the right hand side vectors B for the system of */
/* > linear equations. */
/* > On exit, the solution vectors, X. */
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
/* > \ingroup realPTcomputational */
/* ===================================================================== */
/* Subroutine */
void sptts2_(integer *n, integer *nrhs, real *d__, real *e, real *b, integer *ldb)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sptts2 inputs: n %" FLA_IS ", nrhs %" FLA_IS ", ldb %" FLA_IS "", *n,
             *nrhs, *ldb);
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2;
    real r__1;
    /* Local variables */
    integer i__, j;
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    --d__;
    --e;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    if(*n <= 1)
    {
        if(*n == 1)
        {
            r__1 = 1.f / d__[1];
            sscal_(nrhs, &r__1, &b[b_offset], ldb);
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Solve A * X = B using the factorization A = L*D*L**T, */
    /* overwriting each right hand side vector with its solution. */
    i__1 = *nrhs;
    for(j = 1; j <= i__1; ++j)
    {
        /* Solve L * x = b. */
        i__2 = *n;
        for(i__ = 2; i__ <= i__2; ++i__)
        {
            b[i__ + j * b_dim1] -= b[i__ - 1 + j * b_dim1] * e[i__ - 1];
            /* L10: */
        }
        /* Solve D * L**T * x = b. */
        b[*n + j * b_dim1] /= d__[*n];
        for(i__ = *n - 1; i__ >= 1; --i__)
        {
            b[i__ + j * b_dim1] = b[i__ + j * b_dim1] / d__[i__] - b[i__ + 1 + j * b_dim1] * e[i__];
            /* L20: */
        }
        /* L30: */
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SPTTS2 */
}
/* sptts2_ */
