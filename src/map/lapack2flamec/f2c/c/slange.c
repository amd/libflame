/* slange.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
/*
 *     Modifications Copyright (c) 2025 Advanced Micro Devices, Inc.  All rights reserved.
 */
#include "FLA_f2c.h" /* Table of constant values */
#include "fla_lapack_x86_common.h"

static integer c__1 = 1;
/* > \brief \b SLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest
 * absolute value of any element of a general rectangular matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLANGE + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slange.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slange.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slange.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION SLANGE( NORM, M, N, A, LDA, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER NORM */
/* INTEGER LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANGE returns the value of the one norm, or the Frobenius norm, or */
/* > the infinity norm, or the element of largest absolute value of a */
/* > real matrix A. */
/* > \endverbatim */
/* > */
/* > \return SLANGE */
/* > \verbatim */
/* > */
/* > SLANGE = ( fla_max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* > Specifies the value to be returned in SLANGE as described */
/* > above. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. When M = 0, */
/* > SLANGE is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0. When N = 0, */
/* > SLANGE is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > The m by n matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(M,1). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)), */
/* > where LWORK >= M when NORM = 'I';
otherwise, WORK is not */
/* > referenced. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup realGEauxiliary */
/* ===================================================================== */
real slange_(char *norm, integer *m, integer *n, real *a, integer *lda, real *work)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slange inputs: norm %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "",
                      *norm, *m, *n, *lda);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real ret_val, r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, j_a_dim;
    real sum, temp, scale;
    extern logical lsame_(char *, char *, integer, integer);
    real value;
    extern /* Subroutine */
        void
        slassq_(integer *, real *, integer *, real *, real *);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    /* Function Body */
    value = 0.f;

    /* initialize AOCL context */
    aocl_fla_init();

    if(fla_min(*m, *n) == 0)
    {
        value = 0.f;
    }
    else if(lsame_(norm, "M", 1, 1))
    {
        /* Find fla_max(abs(A(i,j))). */
        value = 0.f;
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *m;
            j_a_dim = j * a_dim1;

#if FLA_ENABLE_AMD_OPT
            /* Select optimized path for AMD architecture*/
            temp = fla_get_max_sabs_element_vector(i__2, a, j_a_dim);

            if(value < temp || temp != temp)
                value = temp;
#else
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                temp = (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1));
                if(value < temp || temp != temp)
                {
                    value = temp;
                }
                /* L10: */
            }
#endif
            /* L20: */
        }
    }
    else if(lsame_(norm, "O", 1, 1) || *(unsigned char *)norm == '1')
    {
        /* Find norm1(A). */
        value = 0.f;
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            sum = 0.f;
            i__2 = *m;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                sum += (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1));
                /* L30: */
            }
            if(value < sum || sum != sum)
            {
                value = sum;
            }
            /* L40: */
        }
    }
    else if(lsame_(norm, "I", 1, 1))
    {
        /* Find normI(A). */
        i__1 = *m;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            work[i__] = 0.f;
            /* L50: */
        }
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *m;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                work[i__] += (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1));
                /* L60: */
            }
            /* L70: */
        }
        value = 0.f;
        i__1 = *m;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            temp = work[i__];
            if(value < temp || temp != temp)
            {
                value = temp;
            }
            /* L80: */
        }
    }
    else if(lsame_(norm, "F", 1, 1) || lsame_(norm, "E", 1, 1))
    {
        /* Find normF(A). */
        scale = 0.f;
        sum = 1.f;
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            slassq_(m, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
            /* L90: */
        }
        value = scale * sqrt(sum);
    }
    ret_val = value;
    AOCL_DTL_TRACE_LOG_EXIT
    return ret_val;
    /* End of SLANGE */
}
/* slange_ */
