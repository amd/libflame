/* ../netlib/spstf2.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static real c_b16 = -1.f;
static real c_b18 = 1.f;
/* > \brief \b SPSTF2 computes the Cholesky factorization with complete pivoting of a real symmetric
 * positive semidefinite matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SPSTF2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spstf2.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spstf2.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spstf2.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SPSTF2( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* REAL TOL */
/* INTEGER INFO, LDA, N, RANK */
/* CHARACTER UPLO */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), WORK( 2*N ) */
/* INTEGER PIV( N ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPSTF2 computes the Cholesky factorization with complete */
/* > pivoting of a real symmetric positive semidefinite matrix A. */
/* > */
/* > The factorization has the form */
/* > P**T * A * P = U**T * U , if UPLO = 'U', */
/* > P**T * A * P = L * L**T, if UPLO = 'L', */
/* > where U is an upper triangular matrix and L is lower triangular, and */
/* > P is stored as vector PIV. */
/* > */
/* > This algorithm does not attempt to check that A is positive */
/* > semidefinite. This version of the algorithm calls level 2 BLAS. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the upper or lower triangular part of the */
/* > symmetric matrix A is stored. */
/* > = 'U': Upper triangular */
/* > = 'L': Lower triangular */
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
/* > n by n upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading n by n lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > */
/* > On exit, if INFO = 0, the factor U or L from the Cholesky */
/* > factorization as above. */
/* > \endverbatim */
/* > */
/* > \param[out] PIV */
/* > \verbatim */
/* > PIV is INTEGER array, dimension (N) */
/* > PIV is such that the nonzero entries are P( PIV(K), K ) = 1. */
/* > \endverbatim */
/* > */
/* > \param[out] RANK */
/* > \verbatim */
/* > RANK is INTEGER */
/* > The rank of A given by the number of steps the algorithm */
/* > completed. */
/* > \endverbatim */
/* > */
/* > \param[in] TOL */
/* > \verbatim */
/* > TOL is REAL */
/* > User defined tolerance. If TOL < 0, then N*U*MAX( A( K,K ) ) */
/* > will be used. The algorithm terminates at the (K-1)st step */
/* > if the pivot <= TOL. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (2*N) */
/* > Work space. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > < 0: If INFO = -K, the K-th argument had an illegal value, */
/* > = 0: algorithm completed successfully, and */
/* > > 0: the matrix A is either rank deficient with computed rank */
/* > as returned in RANK, or is not positive semidefinite. See */
/* > Section 7 of LAPACK Working Note #161 for further */
/* > information. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup realOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
void spstf2_(char *uplo, integer *n, real *a, integer *lda, integer *piv, integer *rank, real *tol,
             real *work, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("spstf2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n,
             *lda);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j;
    real ajj;
    integer pvt;
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *);
    integer itemp;
    extern /* Subroutine */
        void
        sgemv_(char *, integer *, integer *, real *, real *, integer *, real *, integer *, real *,
               real *, integer *);
    real stemp;
    logical upper;
    extern /* Subroutine */
        void
        sswap_(integer *, real *, integer *, real *, integer *);
    real sstop;
    extern real slamch_(char *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer smaxloc_(real *, integer *);
    extern logical sisnan_(real *);
    /* -- LAPACK computational routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* December 2016 */
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
    --work;
    --piv;
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
        xerbla_("SPSTF2", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Initialize PIV */
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        piv[i__] = i__;
        /* L100: */
    }
    /* Compute stopping value */
    pvt = 1;
    ajj = a[pvt + pvt * a_dim1];
    i__1 = *n;
    for(i__ = 2; i__ <= i__1; ++i__)
    {
        if(a[i__ + i__ * a_dim1] > ajj)
        {
            pvt = i__;
            ajj = a[pvt + pvt * a_dim1];
        }
    }
    if(ajj <= 0.f || sisnan_(&ajj))
    {
        *rank = 0;
        *info = 1;
        goto L170;
    }
    /* Compute stopping value if not supplied */
    if(*tol < 0.f)
    {
        sstop = *n * slamch_("Epsilon") * ajj;
    }
    else
    {
        sstop = *tol;
    }
    /* Set first half of WORK to zero, holds dot products */
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        work[i__] = 0.f;
        /* L110: */
    }
    if(upper)
    {
        /* Compute the Cholesky factorization P**T * A * P = U**T * U */
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            /* Find pivot, test for exit, else swap rows and columns */
            /* Update dot products, compute possible pivots which are */
            /* stored in the second half of WORK */
            i__2 = *n;
            for(i__ = j; i__ <= i__2; ++i__)
            {
                if(j > 1)
                {
                    /* Computing 2nd power */
                    r__1 = a[j - 1 + i__ * a_dim1];
                    work[i__] += r__1 * r__1;
                }
                work[*n + i__] = a[i__ + i__ * a_dim1] - work[i__];
                /* L120: */
            }
            if(j > 1)
            {
                i__2 = *n - j + 1;
                itemp = smaxloc_(&work[*n + j], &i__2);
                pvt = itemp + j - 1;
                ajj = work[*n + pvt];
                if(ajj <= sstop || sisnan_(&ajj))
                {
                    a[j + j * a_dim1] = ajj;
                    goto L160;
                }
            }
            if(j != pvt)
            {
                /* Pivot OK, so can now swap pivot rows and columns */
                a[pvt + pvt * a_dim1] = a[j + j * a_dim1];
                i__2 = j - 1;
                sswap_(&i__2, &a[j * a_dim1 + 1], &c__1, &a[pvt * a_dim1 + 1], &c__1);
                if(pvt < *n)
                {
                    i__2 = *n - pvt;
                    sswap_(&i__2, &a[j + (pvt + 1) * a_dim1], lda, &a[pvt + (pvt + 1) * a_dim1],
                           lda);
                }
                i__2 = pvt - j - 1;
                sswap_(&i__2, &a[j + (j + 1) * a_dim1], lda, &a[j + 1 + pvt * a_dim1], &c__1);
                /* Swap dot products and PIV */
                stemp = work[j];
                work[j] = work[pvt];
                work[pvt] = stemp;
                itemp = piv[pvt];
                piv[pvt] = piv[j];
                piv[j] = itemp;
            }
            ajj = sqrt(ajj);
            a[j + j * a_dim1] = ajj;
            /* Compute elements J+1:N of row J */
            if(j < *n)
            {
                i__2 = j - 1;
                i__3 = *n - j;
                sgemv_("Trans", &i__2, &i__3, &c_b16, &a[(j + 1) * a_dim1 + 1], lda,
                       &a[j * a_dim1 + 1], &c__1, &c_b18, &a[j + (j + 1) * a_dim1], lda);
                i__2 = *n - j;
                r__1 = 1.f / ajj;
                sscal_(&i__2, &r__1, &a[j + (j + 1) * a_dim1], lda);
            }
            /* L130: */
        }
    }
    else
    {
        /* Compute the Cholesky factorization P**T * A * P = L * L**T */
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            /* Find pivot, test for exit, else swap rows and columns */
            /* Update dot products, compute possible pivots which are */
            /* stored in the second half of WORK */
            i__2 = *n;
            for(i__ = j; i__ <= i__2; ++i__)
            {
                if(j > 1)
                {
                    /* Computing 2nd power */
                    r__1 = a[i__ + (j - 1) * a_dim1];
                    work[i__] += r__1 * r__1;
                }
                work[*n + i__] = a[i__ + i__ * a_dim1] - work[i__];
                /* L140: */
            }
            if(j > 1)
            {
                i__2 = *n - j + 1;
                itemp = smaxloc_(&work[*n + j], &i__2);
                pvt = itemp + j - 1;
                ajj = work[*n + pvt];
                if(ajj <= sstop || sisnan_(&ajj))
                {
                    a[j + j * a_dim1] = ajj;
                    goto L160;
                }
            }
            if(j != pvt)
            {
                /* Pivot OK, so can now swap pivot rows and columns */
                a[pvt + pvt * a_dim1] = a[j + j * a_dim1];
                i__2 = j - 1;
                sswap_(&i__2, &a[j + a_dim1], lda, &a[pvt + a_dim1], lda);
                if(pvt < *n)
                {
                    i__2 = *n - pvt;
                    sswap_(&i__2, &a[pvt + 1 + j * a_dim1], &c__1, &a[pvt + 1 + pvt * a_dim1],
                           &c__1);
                }
                i__2 = pvt - j - 1;
                sswap_(&i__2, &a[j + 1 + j * a_dim1], &c__1, &a[pvt + (j + 1) * a_dim1], lda);
                /* Swap dot products and PIV */
                stemp = work[j];
                work[j] = work[pvt];
                work[pvt] = stemp;
                itemp = piv[pvt];
                piv[pvt] = piv[j];
                piv[j] = itemp;
            }
            ajj = sqrt(ajj);
            a[j + j * a_dim1] = ajj;
            /* Compute elements J+1:N of column J */
            if(j < *n)
            {
                i__2 = *n - j;
                i__3 = j - 1;
                sgemv_("No Trans", &i__2, &i__3, &c_b16, &a[j + 1 + a_dim1], lda, &a[j + a_dim1],
                       lda, &c_b18, &a[j + 1 + j * a_dim1], &c__1);
                i__2 = *n - j;
                r__1 = 1.f / ajj;
                sscal_(&i__2, &r__1, &a[j + 1 + j * a_dim1], &c__1);
            }
            /* L150: */
        }
    }
    /* Ran to completion, A has full rank */
    *rank = *n;
    goto L170;
L160: /* Rank is number of steps completed. Set INFO = 1 to signal */
    /* that the factorization cannot be used to solve a system. */
    *rank = j - 1;
    *info = 1;
L170:
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SPSTF2 */
}
/* spstf2_ */
