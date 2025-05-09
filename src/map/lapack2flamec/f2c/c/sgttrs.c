/* ../netlib/sgttrs.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b SGTTRS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGTTRS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgttrs.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgttrs.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgttrs.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER INFO, LDB, N, NRHS */
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
/* > SGTTRS solves one of the systems of equations */
/* > A*X = B or A**T*X = B, */
/* > with a tridiagonal matrix A using the LU factorization computed */
/* > by SGTTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the form of the system of equations. */
/* > = 'N': A * X = B (No transpose) */
/* > = 'T': A**T* X = B (Transpose) */
/* > = 'C': A**T* X = B (Conjugate transpose = Transpose) */
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
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
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
void sgttrs_(char *trans, integer *n, integer *nrhs, real *dl, real *d__, real *du, real *du2,
             integer *ipiv, real *b, integer *ldb, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgttrs inputs: trans %c, n %" FLA_IS ", nrhs %" FLA_IS ", ldb %" FLA_IS "",
                      *trans, *n, *nrhs, *ldb);
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3;
    /* Local variables */
    integer j, jb, nb;
    extern /* Subroutine */
        void
        sgtts2_(integer *, integer *, integer *, real *, real *, real *, real *, integer *, real *,
                integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer itrans;
    logical notran;
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
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
    *info = 0;
    notran = *(unsigned char *)trans == 'N' || *(unsigned char *)trans == 'n';
    if(!notran && !(*(unsigned char *)trans == 'T' || *(unsigned char *)trans == 't')
       && !(*(unsigned char *)trans == 'C' || *(unsigned char *)trans == 'c'))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*nrhs < 0)
    {
        *info = -3;
    }
    else if(*ldb < fla_max(*n, 1))
    {
        *info = -10;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGTTRS", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0 || *nrhs == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Decode TRANS */
    if(notran)
    {
        itrans = 0;
    }
    else
    {
        itrans = 1;
    }
    /* Determine the number of right-hand sides to solve at a time. */
    if(*nrhs == 1)
    {
        nb = 1;
    }
    else
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = ilaenv_(&c__1, "SGTTRS", trans, n, nrhs, &c_n1, &c_n1); // , expr subst
        nb = fla_max(i__1, i__2);
    }
    if(nb >= *nrhs)
    {
        sgtts2_(&itrans, n, nrhs, &dl[1], &d__[1], &du[1], &du2[1], &ipiv[1], &b[b_offset], ldb);
    }
    else
    {
        i__1 = *nrhs;
        i__2 = nb;
        for(j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
        {
            /* Computing MIN */
            i__3 = *nrhs - j + 1;
            jb = fla_min(i__3, nb);
            sgtts2_(&itrans, n, &jb, &dl[1], &d__[1], &du[1], &du2[1], &ipiv[1], &b[j * b_dim1 + 1],
                    ldb);
            /* L10: */
        }
    }
    /* End of SGTTRS */
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
/* sgttrs_ */
