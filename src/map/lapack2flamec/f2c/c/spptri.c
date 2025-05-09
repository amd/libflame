/* ../netlib/spptri.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b8 = 1.f;
static integer c__1 = 1;
/* > \brief \b SPPTRI */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SPPTRI + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spptri.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spptri.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spptri.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SPPTRI( UPLO, N, AP, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, N */
/* .. */
/* .. Array Arguments .. */
/* REAL AP( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPPTRI computes the inverse of a real symmetric positive definite */
/* > matrix A using the Cholesky factorization A = U**T*U or A = L*L**T */
/* > computed by SPPTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangular factor is stored in AP;
 */
/* > = 'L': Lower triangular factor is stored in AP. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* > AP is REAL array, dimension (N*(N+1)/2) */
/* > On entry, the triangular factor U or L from the Cholesky */
/* > factorization A = U**T*U or A = L*L**T, packed columnwise as */
/* > a linear array. The j-th column of U or L is stored in the */
/* > array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n. */
/* > */
/* > On exit, the upper or lower triangle of the (symmetric) */
/* > inverse of A, overwriting the input factor U or L. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, the (i,i) element of the factor U or L is */
/* > zero, and the inverse could not be computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
void spptri_(char *uplo, integer *n, real *ap, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("spptri inputs: uplo %c, n %" FLA_IS "", *uplo, *n);
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
    integer j, jc, jj;
    real ajj;
    integer jjn;
    extern real sdot_(integer *, real *, integer *, real *, integer *);
    extern /* Subroutine */
        void
        sspr_(char *, integer *, real *, real *, integer *, real *);
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *);
    logical upper;
    extern /* Subroutine */
        void
        stpmv_(char *, char *, char *, integer *, real *, real *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len),
        stptri_(char *, char *, integer *, real *, integer *);
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
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
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --ap;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SPPTRI", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Invert the triangular Cholesky factor U or L. */
    stptri_(uplo, "Non-unit", n, &ap[1], info);
    if(*info > 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(upper)
    {
        /* Compute the product inv(U) * inv(U)**T. */
        jj = 0;
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            jc = jj + 1;
            jj += j;
            if(j > 1)
            {
                i__2 = j - 1;
                sspr_("Upper", &i__2, &c_b8, &ap[jc], &c__1, &ap[1]);
            }
            ajj = ap[jj];
            sscal_(&j, &ajj, &ap[jc], &c__1);
            /* L10: */
        }
    }
    else
    {
        /* Compute the product inv(L)**T * inv(L). */
        jj = 1;
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            jjn = jj + *n - j + 1;
            i__2 = *n - j + 1;
            ap[jj] = sdot_(&i__2, &ap[jj], &c__1, &ap[jj], &c__1);
            if(j < *n)
            {
                i__2 = *n - j;
                stpmv_("Lower", "Transpose", "Non-unit", &i__2, &ap[jjn], &ap[jj + 1], &c__1);
            }
            jj = jjn;
            /* L20: */
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SPPTRI */
}
/* spptri_ */
