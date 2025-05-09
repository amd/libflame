/* ../netlib/clacp2.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLACP2 copies all or part of a real two-dimensional array to a complex array. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLACP2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacp2.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacp2.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacp2.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLACP2( UPLO, M, N, A, LDA, B, LDB ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER LDA, LDB, M, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ) */
/* COMPLEX B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLACP2 copies all or part of a real two-dimensional matrix A to a */
/* > complex matrix B. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies the part of the matrix A to be copied to B. */
/* > = 'U': Upper triangular part */
/* > = 'L': Lower triangular part */
/* > Otherwise: All of the matrix A */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > The m by n matrix A. If UPLO = 'U', only the upper trapezium */
/* > is accessed;
if UPLO = 'L', only the lower trapezium is */
/* > accessed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,N) */
/* > On exit, B = A in the locations specified by UPLO. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,M). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
void clacp2_(char *uplo, integer *m, integer *n, real *a, integer *lda, complex *b, integer *ldb)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clacp2 inputs: uplo %c, m %lld, n %lld, lda %lld, ldb %lld", *uplo, *m,
             *n, *lda, *ldb);
#else
    snprintf(buffer, 256, "clacp2 inputs: uplo %c, m %d, n %d, lda %d, ldb %d", *uplo, *m, *n, *lda,
             *ldb);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    integer i__, j;
    extern logical lsame_(char *, char *, integer, integer);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    if(lsame_(uplo, "U", 1, 1))
    {
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = fla_min(j, *m);
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                i__3 = i__ + j * b_dim1;
                i__4 = i__ + j * a_dim1;
                b[i__3].r = a[i__4];
                b[i__3].i = 0.f; // , expr subst
                /* L10: */
            }
            /* L20: */
        }
    }
    else if(lsame_(uplo, "L", 1, 1))
    {
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *m;
            for(i__ = j; i__ <= i__2; ++i__)
            {
                i__3 = i__ + j * b_dim1;
                i__4 = i__ + j * a_dim1;
                b[i__3].r = a[i__4];
                b[i__3].i = 0.f; // , expr subst
                /* L30: */
            }
            /* L40: */
        }
    }
    else
    {
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *m;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                i__3 = i__ + j * b_dim1;
                i__4 = i__ + j * a_dim1;
                b[i__3].r = a[i__4];
                b[i__3].i = 0.f; // , expr subst
                /* L50: */
            }
            /* L60: */
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLACP2 */
}
/* clacp2_ */
