/* ../netlib/cptcon.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CPTCON */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CPTCON + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cptcon.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cptcon.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cptcon.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CPTCON( N, D, E, ANORM, RCOND, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, N */
/* REAL ANORM, RCOND */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ), RWORK( * ) */
/* COMPLEX E( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPTCON computes the reciprocal of the condition number (in the */
/* > 1-norm) of a complex Hermitian positive definite tridiagonal matrix */
/* > using the factorization A = L*D*L**H or A = U**H*D*U computed by */
/* > CPTTRF. */
/* > */
/* > Norm(inv(A)) is computed by a direct method, and the reciprocal of */
/* > the condition number is computed as */
/* > RCOND = 1 / (ANORM * norm(inv(A))). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The n diagonal elements of the diagonal matrix D from the */
/* > factorization of A, as computed by CPTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is COMPLEX array, dimension (N-1) */
/* > The (n-1) off-diagonal elements of the unit bidiagonal factor */
/* > U or L from the factorization of A, as computed by CPTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* > ANORM is REAL */
/* > The 1-norm of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* > RCOND is REAL */
/* > The reciprocal of the condition number of the matrix A, */
/* > computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is the */
/* > 1-norm of inv(A) computed in this routine. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (N) */
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
/* > \ingroup complexPTcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The method used is described in Nicholas J. Higham, "Efficient */
/* > Algorithms for Computing the Condition Number of a Tridiagonal */
/* > Matrix", SIAM J. Sci. Stat. Comput., Vol. 7, No. 1, January 1986. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void cptcon_(integer *n, real *d__, complex *e, real *anorm, real *rcond, real *rwork,
             integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "cptcon inputs: n %lld", *n);
#else
    snprintf(buffer, 256, "cptcon inputs: n %d", *n);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer i__1;
    real r__1;
    /* Builtin functions */
    double c_abs(complex *);
    /* Local variables */
    integer i__, ix;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer isamax_(integer *, real *, integer *);
    real ainvnm;
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments. */
    /* Parameter adjustments */
    --rwork;
    --e;
    --d__;
    /* Function Body */
    *info = 0;
    if(*n < 0)
    {
        *info = -1;
    }
    else if(*anorm < 0.f)
    {
        *info = -4;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CPTCON", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    *rcond = 0.f;
    if(*n == 0)
    {
        *rcond = 1.f;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    else if(*anorm == 0.f)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Check that D(1:N) is positive. */
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if(d__[i__] <= 0.f)
        {
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
        /* L10: */
    }
    /* Solve M(A) * x = e, where M(A) = (m(i,j)) is given by */
    /* m(i,j) = f2c_abs(A(i,j)); i = j; */
    /* m(i,j) = -f2c_abs(A(i,j)), i .ne. j, */
    /* and e = [ 1, 1, ..., 1 ]**T. Note M(A) = M(L)*D*M(L)**H. */
    /* Solve M(L) * x = e. */
    rwork[1] = 1.f;
    i__1 = *n;
    for(i__ = 2; i__ <= i__1; ++i__)
    {
        rwork[i__] = rwork[i__ - 1] * c_abs(&e[i__ - 1]) + 1.f;
        /* L20: */
    }
    /* Solve D * M(L)**H * x = b. */
    rwork[*n] /= d__[*n];
    for(i__ = *n - 1; i__ >= 1; --i__)
    {
        rwork[i__] = rwork[i__] / d__[i__] + rwork[i__ + 1] * c_abs(&e[i__]);
        /* L30: */
    }
    /* Compute AINVNM = fla_max(x(i)), 1<=i<=n. */
    ix = isamax_(n, &rwork[1], &c__1);
    ainvnm = (r__1 = rwork[ix], f2c_abs(r__1));
    /* Compute the reciprocal condition number. */
    if(ainvnm != 0.f)
    {
        *rcond = 1.f / ainvnm / *anorm;
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CPTCON */
}
/* cptcon_ */
