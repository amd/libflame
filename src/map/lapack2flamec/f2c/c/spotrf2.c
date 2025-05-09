/* ../netlib/v3.9.0/spotrf2.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b9 = 1.f;
static real c_b11 = -1.f;
/* > \brief \b SPOTRF2 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* Definition: */
/* =========== */
/* SUBROUTINE SPOTRF2( UPLO, N, A, LDA, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPOTRF2 computes the Cholesky factorization of a real symmetric */
/* > positive definite matrix A using the recursive algorithm. */
/* > */
/* > The factorization has the form */
/* > A = U**T * U, if UPLO = 'U', or */
/* > A = L * L**T, if UPLO = 'L', */
/* > where U is an upper triangular matrix and L is lower triangular. */
/* > */
/* > This is the recursive version of the algorithm. It divides */
/* > the matrix into four submatrices: */
/* > */
/* > [ A11 | A12 ] where A11 is n1 by n1 and A22 is n2 by n2 */
/* > A = [ -----|----- ] with n1 = n/2 */
/* > [ A21 | A22 ] n2 = n-n1 */
/* > */
/* > The subroutine calls itself to factor A11. Update and scale A21 */
/* > or A12, update A22 then call itself to factor A22. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
 */
/* > = 'L': Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the symmetric matrix A. If UPLO = 'U', the leading */
/* > N-by-N upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading N-by-N lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > */
/* > On exit, if INFO = 0, the factor U or L from the Cholesky */
/* > factorization A = U**T*U or A = L*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, the leading minor of order i is not */
/* > positive definite, and the factorization could not be */
/* > completed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2017 */
/* > \ingroup realPOcomputational */
/* ===================================================================== */
/* Subroutine */
void spotrf2_(char *uplo, integer *n, real *a, integer *lda, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer n1, n2;
    extern logical lsame_(char *, char *, integer, integer);
    integer iinfo;
    logical upper;
    extern /* Subroutine */
        void
        strsm_(char *, char *, char *, char *, integer *, integer *, real *, real *, integer *,
               real *, integer *),
        ssyrk_(char *, char *, integer *, integer *, real *, real *, integer *, real *, real *,
               integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern logical sisnan_(real *);
    /* -- LAPACK computational routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2017 */
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
    /* Test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
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
    else if(*lda < fla_max(1, *n))
    {
        *info = -4;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SPOTRF2", &i__1, (ftnlen)7);
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        return;
    }
    /* N=1 case */
    if(*n == 1)
    {
        /* Test for non-positive-definiteness */
        if(a[a_dim1 + 1] <= 0.f || sisnan_(&a[a_dim1 + 1]))
        {
            *info = 1;
            return;
        }
        /* Factor */
        a[a_dim1 + 1] = sqrt(a[a_dim1 + 1]);
        /* Use recursive code */
    }
    else
    {
        n1 = *n / 2;
        n2 = *n - n1;
        /* Factor A11 */
        spotrf2_(uplo, &n1, &a[a_dim1 + 1], lda, &iinfo);
        if(iinfo != 0)
        {
            *info = iinfo;
            return;
        }
        /* Compute the Cholesky factorization A = U**T*U */
        if(upper)
        {
            /* Update and scale A12 */
            strsm_("L", "U", "T", "N", &n1, &n2, &c_b9, &a[a_dim1 + 1], lda,
                   &a[(n1 + 1) * a_dim1 + 1], lda);
            /* Update and factor A22 */
            ssyrk_(uplo, "T", &n2, &n1, &c_b11, &a[(n1 + 1) * a_dim1 + 1], lda, &c_b9,
                   &a[n1 + 1 + (n1 + 1) * a_dim1], lda);
            spotrf2_(uplo, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], lda, &iinfo);
            if(iinfo != 0)
            {
                *info = iinfo + n1;
                return;
            }
            /* Compute the Cholesky factorization A = L*L**T */
        }
        else
        {
            /* Update and scale A21 */
            strsm_("R", "L", "T", "N", &n2, &n1, &c_b9, &a[a_dim1 + 1], lda, &a[n1 + 1 + a_dim1],
                   lda);
            /* Update and factor A22 */
            ssyrk_(uplo, "N", &n2, &n1, &c_b11, &a[n1 + 1 + a_dim1], lda, &c_b9,
                   &a[n1 + 1 + (n1 + 1) * a_dim1], lda);
            spotrf2_(uplo, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], lda, &iinfo);
            if(iinfo != 0)
            {
                *info = iinfo + n1;
                return;
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SPOTRF2 */
}
/* spotrf2_ */
