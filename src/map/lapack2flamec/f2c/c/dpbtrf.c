/* ../netlib/dpbtrf.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b18 = 1.;
static doublereal c_b21 = -1.;
static integer c__33 = 33;
/* > \brief \b DPBTRF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DPBTRF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpbtrf.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpbtrf.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpbtrf.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, KD, LDAB, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION AB( LDAB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPBTRF computes the Cholesky factorization of a real symmetric */
/* > positive definite band matrix A. */
/* > */
/* > The factorization has the form */
/* > A = U**T * U, if UPLO = 'U', or */
/* > A = L * L**T, if UPLO = 'L', */
/* > where U is an upper triangular matrix and L is lower triangular. */
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
/* > \param[in] KD */
/* > \verbatim */
/* > KD is INTEGER */
/* > The number of superdiagonals of the matrix A if UPLO = 'U', */
/* > or the number of subdiagonals if UPLO = 'L'. KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* > AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* > On entry, the upper or lower triangle of the symmetric band */
/* > matrix A, stored in the first KD+1 rows of the array. The */
/* > j-th column of A is stored in the j-th column of the array AB */
/* > as follows: */
/* > if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j;
 */
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=fla_min(n,j+kd). */
/* > */
/* > On exit, if INFO = 0, the triangular factor U or L from the */
/* > Cholesky factorization A = U**T*U or A = L*L**T of the band */
/* > matrix A, in the same storage format as A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KD+1. */
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
/* > \date November 2011 */
/* > \ingroup doubleOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The band storage scheme is illustrated by the following example, when */
/* > N = 6, KD = 2, and UPLO = 'U': */
/* > */
/* > On entry: On exit: */
/* > */
/* > * * a13 a24 a35 a46 * * u13 u24 u35 u46 */
/* > * a12 a23 a34 a45 a56 * u12 u23 u34 u45 u56 */
/* > a11 a22 a33 a44 a55 a66 u11 u22 u33 u44 u55 u66 */
/* > */
/* > Similarly, if UPLO = 'L' the format of A is as follows: */
/* > */
/* > On entry: On exit: */
/* > */
/* > a11 a22 a33 a44 a55 a66 l11 l22 l33 l44 l55 l66 */
/* > a21 a32 a43 a54 a65 * l21 l32 l43 l54 l65 * */
/* > a31 a42 a53 a64 * * l31 l42 l53 l64 * * */
/* > */
/* > Array elements marked * are not used by the routine. */
/* > \endverbatim */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989 */
/* ===================================================================== */
/* Subroutine */
void dpbtrf_(char *uplo, integer *n, integer *kd, doublereal *ab, integer *ldab, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dpbtrf inputs: uplo %c, n %" FLA_IS ", kd %" FLA_IS ", ldab %" FLA_IS "",
                      *uplo, *n, *kd, *ldab);
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    integer i__, j, i2, i3, ib, nb, ii, jj;
    doublereal work[1056] /* was [33][32] */
        ;
    extern /* Subroutine */
        void
        dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *),
        dsyrk_(char *, char *, integer *, integer *, doublereal *, doublereal *, integer *,
               doublereal *, doublereal *, integer *),
        dpbtf2_(char *, integer *, integer *, doublereal *, integer *, integer *),
        dpotf2_(char *, integer *, doublereal *, integer *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
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
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    /* Function Body */
    *info = 0;
    if(!lsame_(uplo, "U", 1, 1) && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*kd < 0)
    {
        *info = -3;
    }
    else if(*ldab < *kd + 1)
    {
        *info = -5;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DPBTRF", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Determine the block size for this environment */
    nb = ilaenv_(&c__1, "DPBTRF", uplo, n, kd, &c_n1, &c_n1);
    /* The block size must not exceed the semi-bandwidth KD, and must not */
    /* exceed the limit set by the size of the local array WORK. */
    nb = fla_min(nb, 32);
    if(nb <= 1 || nb > *kd)
    {
        /* Use unblocked code */
        dpbtf2_(uplo, n, kd, &ab[ab_offset], ldab, info);
    }
    else
    {
        /* Use blocked code */
        if(lsame_(uplo, "U", 1, 1))
        {
            /* Compute the Cholesky factorization of a symmetric band */
            /* matrix, given the upper triangle of the matrix in band */
            /* storage. */
            /* Zero the upper triangle of the work array. */
            i__1 = nb;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = j - 1;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    work[i__ + j * 33 - 34] = 0.;
                    /* L10: */
                }
                /* L20: */
            }
            /* Process the band matrix one diagonal block at a time. */
            i__1 = *n;
            i__2 = nb;
            for(i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
            {
                /* Computing MIN */
                i__3 = nb;
                i__4 = *n - i__ + 1; // , expr subst
                ib = fla_min(i__3, i__4);
                /* Factorize the diagonal block */
                i__3 = *ldab - 1;
                dpotf2_(uplo, &ib, &ab[*kd + 1 + i__ * ab_dim1], &i__3, &ii);
                if(ii != 0)
                {
                    *info = i__ + ii - 1;
                    goto L150;
                }
                if(i__ + ib <= *n)
                {
                    /* Update the relevant part of the trailing submatrix. */
                    /* If A11 denotes the diagonal block which has just been */
                    /* factorized, then we need to update the remaining */
                    /* blocks in the diagram: */
                    /* A11 A12 A13 */
                    /* A22 A23 */
                    /* A33 */
                    /* The numbers of rows and columns in the partitioning */
                    /* are IB, I2, I3 respectively. The blocks A12, A22 and */
                    /* A23 are empty if IB = KD. The upper triangle of A13 */
                    /* lies outside the band. */
                    /* Computing MIN */
                    i__3 = *kd - ib;
                    i__4 = *n - i__ - ib + 1; // , expr subst
                    i2 = fla_min(i__3, i__4);
                    /* Computing MIN */
                    i__3 = ib;
                    i__4 = *n - i__ - *kd + 1; // , expr subst
                    i3 = fla_min(i__3, i__4);
                    if(i2 > 0)
                    {
                        /* Update A12 */
                        i__3 = *ldab - 1;
                        i__4 = *ldab - 1;
                        dtrsm_("Left", "Upper", "Transpose", "Non-unit", &ib, &i2, &c_b18,
                               &ab[*kd + 1 + i__ * ab_dim1], &i__3,
                               &ab[*kd + 1 - ib + (i__ + ib) * ab_dim1], &i__4);
                        /* Update A22 */
                        i__3 = *ldab - 1;
                        i__4 = *ldab - 1;
                        dsyrk_("Upper", "Transpose", &i2, &ib, &c_b21,
                               &ab[*kd + 1 - ib + (i__ + ib) * ab_dim1], &i__3, &c_b18,
                               &ab[*kd + 1 + (i__ + ib) * ab_dim1], &i__4);
                    }
                    if(i3 > 0)
                    {
                        /* Copy the lower triangle of A13 into the work array. */
                        i__3 = i3;
                        for(jj = 1; jj <= i__3; ++jj)
                        {
                            i__4 = ib;
                            for(ii = jj; ii <= i__4; ++ii)
                            {
                                work[ii + jj * 33 - 34]
                                    = ab[ii - jj + 1 + (jj + i__ + *kd - 1) * ab_dim1];
                                /* L30: */
                            }
                            /* L40: */
                        }
                        /* Update A13 (in the work array). */
                        i__3 = *ldab - 1;
                        dtrsm_("Left", "Upper", "Transpose", "Non-unit", &ib, &i3, &c_b18,
                               &ab[*kd + 1 + i__ * ab_dim1], &i__3, work, &c__33);
                        /* Update A23 */
                        if(i2 > 0)
                        {
                            i__3 = *ldab - 1;
                            i__4 = *ldab - 1;
                            dgemm_("Transpose", "No Transpose", &i2, &i3, &ib, &c_b21,
                                   &ab[*kd + 1 - ib + (i__ + ib) * ab_dim1], &i__3, work, &c__33,
                                   &c_b18, &ab[ib + 1 + (i__ + *kd) * ab_dim1], &i__4);
                        }
                        /* Update A33 */
                        i__3 = *ldab - 1;
                        dsyrk_("Upper", "Transpose", &i3, &ib, &c_b21, work, &c__33, &c_b18,
                               &ab[*kd + 1 + (i__ + *kd) * ab_dim1], &i__3);
                        /* Copy the lower triangle of A13 back into place. */
                        i__3 = i3;
                        for(jj = 1; jj <= i__3; ++jj)
                        {
                            i__4 = ib;
                            for(ii = jj; ii <= i__4; ++ii)
                            {
                                ab[ii - jj + 1 + (jj + i__ + *kd - 1) * ab_dim1]
                                    = work[ii + jj * 33 - 34];
                                /* L50: */
                            }
                            /* L60: */
                        }
                    }
                }
                /* L70: */
            }
        }
        else
        {
            /* Compute the Cholesky factorization of a symmetric band */
            /* matrix, given the lower triangle of the matrix in band */
            /* storage. */
            /* Zero the lower triangle of the work array. */
            i__2 = nb;
            for(j = 1; j <= i__2; ++j)
            {
                i__1 = nb;
                for(i__ = j + 1; i__ <= i__1; ++i__)
                {
                    work[i__ + j * 33 - 34] = 0.;
                    /* L80: */
                }
                /* L90: */
            }
            /* Process the band matrix one diagonal block at a time. */
            i__2 = *n;
            i__1 = nb;
            for(i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
            {
                /* Computing MIN */
                i__3 = nb;
                i__4 = *n - i__ + 1; // , expr subst
                ib = fla_min(i__3, i__4);
                /* Factorize the diagonal block */
                i__3 = *ldab - 1;
                dpotf2_(uplo, &ib, &ab[i__ * ab_dim1 + 1], &i__3, &ii);
                if(ii != 0)
                {
                    *info = i__ + ii - 1;
                    goto L150;
                }
                if(i__ + ib <= *n)
                {
                    /* Update the relevant part of the trailing submatrix. */
                    /* If A11 denotes the diagonal block which has just been */
                    /* factorized, then we need to update the remaining */
                    /* blocks in the diagram: */
                    /* A11 */
                    /* A21 A22 */
                    /* A31 A32 A33 */
                    /* The numbers of rows and columns in the partitioning */
                    /* are IB, I2, I3 respectively. The blocks A21, A22 and */
                    /* A32 are empty if IB = KD. The lower triangle of A31 */
                    /* lies outside the band. */
                    /* Computing MIN */
                    i__3 = *kd - ib;
                    i__4 = *n - i__ - ib + 1; // , expr subst
                    i2 = fla_min(i__3, i__4);
                    /* Computing MIN */
                    i__3 = ib;
                    i__4 = *n - i__ - *kd + 1; // , expr subst
                    i3 = fla_min(i__3, i__4);
                    if(i2 > 0)
                    {
                        /* Update A21 */
                        i__3 = *ldab - 1;
                        i__4 = *ldab - 1;
                        dtrsm_("Right", "Lower", "Transpose", "Non-unit", &i2, &ib, &c_b18,
                               &ab[i__ * ab_dim1 + 1], &i__3, &ab[ib + 1 + i__ * ab_dim1], &i__4);
                        /* Update A22 */
                        i__3 = *ldab - 1;
                        i__4 = *ldab - 1;
                        dsyrk_("Lower", "No Transpose", &i2, &ib, &c_b21,
                               &ab[ib + 1 + i__ * ab_dim1], &i__3, &c_b18,
                               &ab[(i__ + ib) * ab_dim1 + 1], &i__4);
                    }
                    if(i3 > 0)
                    {
                        /* Copy the upper triangle of A31 into the work array. */
                        i__3 = ib;
                        for(jj = 1; jj <= i__3; ++jj)
                        {
                            i__4 = fla_min(jj, i3);
                            for(ii = 1; ii <= i__4; ++ii)
                            {
                                work[ii + jj * 33 - 34]
                                    = ab[*kd + 1 - jj + ii + (jj + i__ - 1) * ab_dim1];
                                /* L100: */
                            }
                            /* L110: */
                        }
                        /* Update A31 (in the work array). */
                        i__3 = *ldab - 1;
                        dtrsm_("Right", "Lower", "Transpose", "Non-unit", &i3, &ib, &c_b18,
                               &ab[i__ * ab_dim1 + 1], &i__3, work, &c__33);
                        /* Update A32 */
                        if(i2 > 0)
                        {
                            i__3 = *ldab - 1;
                            i__4 = *ldab - 1;
                            dgemm_("No transpose", "Transpose", &i3, &i2, &ib, &c_b21, work, &c__33,
                                   &ab[ib + 1 + i__ * ab_dim1], &i__3, &c_b18,
                                   &ab[*kd + 1 - ib + (i__ + ib) * ab_dim1], &i__4);
                        }
                        /* Update A33 */
                        i__3 = *ldab - 1;
                        dsyrk_("Lower", "No Transpose", &i3, &ib, &c_b21, work, &c__33, &c_b18,
                               &ab[(i__ + *kd) * ab_dim1 + 1], &i__3);
                        /* Copy the upper triangle of A31 back into place. */
                        i__3 = ib;
                        for(jj = 1; jj <= i__3; ++jj)
                        {
                            i__4 = fla_min(jj, i3);
                            for(ii = 1; ii <= i__4; ++ii)
                            {
                                ab[*kd + 1 - jj + ii + (jj + i__ - 1) * ab_dim1]
                                    = work[ii + jj * 33 - 34];
                                /* L120: */
                            }
                            /* L130: */
                        }
                    }
                }
                /* L140: */
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
L150:
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DPBTRF */
}
/* dpbtrf_ */
