/* clanhb.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CLANHB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm,
 * or the ele ment of largest absolute value of a Hermitian band matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLANHB + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhb.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhb.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhb.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION CLANHB( NORM, UPLO, N, K, AB, LDAB, */
/* WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER NORM, UPLO */
/* INTEGER K, LDAB, N */
/* .. */
/* .. Array Arguments .. */
/* REAL WORK( * ) */
/* COMPLEX AB( LDAB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLANHB returns the value of the one norm, or the Frobenius norm, or */
/* > the infinity norm, or the element of largest absolute value of an */
/* > n by n hermitian band matrix A, with k super-diagonals. */
/* > \endverbatim */
/* > */
/* > \return CLANHB */
/* > \verbatim */
/* > */
/* > CLANHB = ( fla_max(abs(A(i,j))), NORM = 'M' or 'm' */
/* > ( */
/* > ( norm1(A), NORM = '1', 'O' or 'o' */
/* > ( */
/* > ( normI(A), NORM = 'I' or 'i' */
/* > ( */
/* > ( normF(A), NORM = 'F', 'f', 'E' or 'e' */
/* > */
/* > where norm1 denotes the one norm of a matrix (maximum column sum), */
/* > normI denotes the infinity norm of a matrix (maximum row sum) and */
/* > normF denotes the Frobenius norm of a matrix (square root of sum of */
/* > squares). Note that fla_max(abs(A(i,j))) is not a consistent matrix norm. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] NORM */
/* > \verbatim */
/* > NORM is CHARACTER*1 */
/* > Specifies the value to be returned in CLANHB as described */
/* > above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the upper or lower triangular part of the */
/* > band matrix A is supplied. */
/* > = 'U': Upper triangular */
/* > = 'L': Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. When N = 0, CLANHB is */
/* > set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of super-diagonals or sub-diagonals of the */
/* > band matrix A. K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* > AB is COMPLEX array, dimension (LDAB,N) */
/* > The upper or lower triangle of the hermitian band matrix A, */
/* > stored in the first K+1 rows of AB. The j-th column of A is */
/* > stored in the j-th column of the array AB as follows: */
/* > if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for fla_max(1,j-k)<=i<=j;
 */
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=fla_min(n,j+k). */
/* > Note that the imaginary parts of the diagonal elements need */
/* > not be set and are assumed to be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= K+1. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)), */
/* > where LWORK >= N when NORM = 'I' or '1' or 'O';
otherwise, */
/* > WORK is not referenced. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
real clanhb_(char *norm, char *uplo, integer *n, integer *k, complex *ab, integer *ldab, real *work)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clanhb inputs: norm %c, uplo %c, n %lld, k %lld, ldab %lld", *norm,
             *uplo, *n, *k, *ldab);
#else
    snprintf(buffer, 256, "clanhb inputs: norm %c, uplo %c, n %d, k %d, ldab %d", *norm, *uplo, *n,
             *k, *ldab);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    real ret_val, r__1;
    /* Builtin functions */
    double c_abs(complex *), sqrt(doublereal);
    /* Local variables */
    integer i__, j, l;
    real sum, absa, scale;
    extern logical lsame_(char *, char *, integer, integer);
    real value;
    extern /* Subroutine */
        void
        classq_(integer *, complex *, integer *, real *, real *);
    extern logical sisnan_(real *);
    /* -- LAPACK auxiliary routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
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
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --work;
    /* Function Body */
    value = 0.f;
    if(*n == 0)
    {
        value = 0.f;
    }
    else if(lsame_(norm, "M", 1, 1))
    {
        /* Find fla_max(abs(A(i,j))). */
        value = 0.f;
        if(lsame_(uplo, "U", 1, 1))
        {
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                /* Computing MAX */
                i__2 = *k + 2 - j;
                i__3 = *k;
                for(i__ = fla_max(i__2, 1); i__ <= i__3; ++i__)
                {
                    sum = c_abs(&ab[i__ + j * ab_dim1]);
                    if(value < sum || sisnan_(&sum))
                    {
                        value = sum;
                    }
                    /* L10: */
                }
                i__3 = *k + 1 + j * ab_dim1;
                sum = (r__1 = ab[i__3].r, f2c_abs(r__1));
                if(value < sum || sisnan_(&sum))
                {
                    value = sum;
                }
                /* L20: */
            }
        }
        else
        {
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                i__3 = j * ab_dim1 + 1;
                sum = (r__1 = ab[i__3].r, f2c_abs(r__1));
                if(value < sum || sisnan_(&sum))
                {
                    value = sum;
                }
                /* Computing MIN */
                i__2 = *n + 1 - j;
                i__4 = *k + 1; // , expr subst
                i__3 = fla_min(i__2, i__4);
                for(i__ = 2; i__ <= i__3; ++i__)
                {
                    sum = c_abs(&ab[i__ + j * ab_dim1]);
                    if(value < sum || sisnan_(&sum))
                    {
                        value = sum;
                    }
                    /* L30: */
                }
                /* L40: */
            }
        }
    }
    else if(lsame_(norm, "I", 1, 1) || lsame_(norm, "O", 1, 1) || *(unsigned char *)norm == '1')
    {
        /* Find normI(A) ( = norm1(A), since A is hermitian). */
        value = 0.f;
        if(lsame_(uplo, "U", 1, 1))
        {
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                sum = 0.f;
                l = *k + 1 - j;
                /* Computing MAX */
                i__3 = 1;
                i__2 = j - *k; // , expr subst
                i__4 = j - 1;
                for(i__ = fla_max(i__3, i__2); i__ <= i__4; ++i__)
                {
                    absa = c_abs(&ab[l + i__ + j * ab_dim1]);
                    sum += absa;
                    work[i__] += absa;
                    /* L50: */
                }
                i__4 = *k + 1 + j * ab_dim1;
                work[j] = sum + (r__1 = ab[i__4].r, f2c_abs(r__1));
                /* L60: */
            }
            i__1 = *n;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                sum = work[i__];
                if(value < sum || sisnan_(&sum))
                {
                    value = sum;
                }
                /* L70: */
            }
        }
        else
        {
            i__1 = *n;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                work[i__] = 0.f;
                /* L80: */
            }
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                i__4 = j * ab_dim1 + 1;
                sum = work[j] + (r__1 = ab[i__4].r, f2c_abs(r__1));
                l = 1 - j;
                /* Computing MIN */
                i__3 = *n;
                i__2 = j + *k; // , expr subst
                i__4 = fla_min(i__3, i__2);
                for(i__ = j + 1; i__ <= i__4; ++i__)
                {
                    absa = c_abs(&ab[l + i__ + j * ab_dim1]);
                    sum += absa;
                    work[i__] += absa;
                    /* L90: */
                }
                if(value < sum || sisnan_(&sum))
                {
                    value = sum;
                }
                /* L100: */
            }
        }
    }
    else if(lsame_(norm, "F", 1, 1) || lsame_(norm, "E", 1, 1))
    {
        /* Find normF(A). */
        scale = 0.f;
        sum = 1.f;
        if(*k > 0)
        {
            if(lsame_(uplo, "U", 1, 1))
            {
                i__1 = *n;
                for(j = 2; j <= i__1; ++j)
                {
                    /* Computing MIN */
                    i__3 = j - 1;
                    i__4 = fla_min(i__3, *k);
                    /* Computing MAX */
                    i__2 = *k + 2 - j;
                    classq_(&i__4, &ab[fla_max(i__2, 1) + j * ab_dim1], &c__1, &scale, &sum);
                    /* L110: */
                }
                l = *k + 1;
            }
            else
            {
                i__1 = *n - 1;
                for(j = 1; j <= i__1; ++j)
                {
                    /* Computing MIN */
                    i__3 = *n - j;
                    i__4 = fla_min(i__3, *k);
                    classq_(&i__4, &ab[j * ab_dim1 + 2], &c__1, &scale, &sum);
                    /* L120: */
                }
                l = 1;
            }
            sum *= 2;
        }
        else
        {
            l = 1;
        }
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__4 = l + j * ab_dim1;
            if(ab[i__4].r != 0.f)
            {
                i__4 = l + j * ab_dim1;
                absa = (r__1 = ab[i__4].r, f2c_abs(r__1));
                if(scale < absa)
                {
                    /* Computing 2nd power */
                    r__1 = scale / absa;
                    sum = sum * (r__1 * r__1) + 1.f;
                    scale = absa;
                }
                else
                {
                    /* Computing 2nd power */
                    r__1 = absa / scale;
                    sum += r__1 * r__1;
                }
            }
            /* L130: */
        }
        value = scale * sqrt(sum);
    }
    ret_val = value;
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return ret_val;
    /* End of CLANHB */
}
/* clanhb_ */
