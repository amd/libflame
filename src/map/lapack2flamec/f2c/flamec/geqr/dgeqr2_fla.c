/* dgeqr2.f -- translated by f2c (version 20000121). You must link the resulting object file with
 * the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DGEQR2 computes the QR factorization of a general rectangular matrix using an
 * unblocked algorit hm. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGEQR2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqr2.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqr2.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqr2.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEQR2 computes a QR factorization of a real m-by-n matrix A: */
/* > */
/* > A = Q * ( R ), */
/* > ( 0 ) */
/* > */
/* > where: */
/* > */
/* > Q is a m-by-m orthogonal matrix;
 */
/* > R is an upper-triangular n-by-n matrix;
 */
/* > 0 is a (m-n)-by-n zero matrix, if m > n. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
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
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the m by n matrix A. */
/* > On exit, the elements on and above the diagonal of the array */
/* > contain the fla_min(m,n) by n upper trapezoidal matrix R (R is */
/* > upper triangular if m >= n);
the elements below the diagonal, */
/* > with the array TAU, represent the orthogonal matrix Q as a */
/* > product of elementary reflectors (see Further Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION array, dimension (fla_min(M,N)) */
/* > The scalar factors of the elementary reflectors (see Further */
/* > Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (N) */
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
/* > \ingroup doubleGEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix Q is represented as a product of elementary reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(k), where k = fla_min(m,n). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**T */
/* > */
/* > where tau is a real scalar, and v is a real vector with */
/* > v(1:i-1) = 0 and v(i) = 1;
v(i+1:m) is stored on exit in A(i+1:m,i), */
/* > and tau in TAU(i). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void dgeqr2_fla(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau,
                doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    integer i__, k;
    extern /* Subroutine */
        void
        dlarf_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
               integer *, doublereal *),
        dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    doublereal aii;
    /* -- LAPACK computational routine -- */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -4;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGEQR2", &i__1, (ftnlen)6);
        return;
    }
    k = fla_min(*m, *n);
    i__1 = k;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        /* Generate elementary reflector H(i) to annihilate A(i+1:m,i) */
        i__2 = *m - i__ + 1;
        /* Computing MIN */
        i__3 = i__ + 1;
        dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[fla_min(i__3, *m) + i__ * a_dim1], &c__1,
                &tau[i__]);
        if(i__ < *n)
        {
            /* Apply H(i) to A(i:m,i+1:n) from the left */
            aii = a[i__ + i__ * a_dim1];
            a[i__ + i__ * a_dim1] = 1.;
            i__2 = *m - i__ + 1;
            i__3 = *n - i__;
            dlarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &tau[i__],
                   &a[i__ + (i__ + 1) * a_dim1], lda, &work[1]);
            a[i__ + i__ * a_dim1] = aii;
        }
        /* L10: */
    }
    return;
    /* End of DGEQR2 */
}
/* dgeqr2_ */
