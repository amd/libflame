/* chptrd.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b2 = {0.f, 0.f};
static integer c__1 = 1;
/* > \brief \b CHPTRD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHPTRD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chptrd.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chptrd.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chptrd.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHPTRD( UPLO, N, AP, D, E, TAU, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, N */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ), E( * ) */
/* COMPLEX AP( * ), TAU( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPTRD reduces a complex Hermitian matrix A stored in packed form to */
/* > real symmetric tridiagonal form T by a unitary similarity */
/* > transformation: Q**H * A * Q = T. */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* > AP is COMPLEX array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the Hermitian matrix */
/* > A, packed columnwise in a linear array. The j-th column of A */
/* > is stored in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > On exit, if UPLO = 'U', the diagonal and first superdiagonal */
/* > of A are overwritten by the corresponding elements of the */
/* > tridiagonal matrix T, and the elements above the first */
/* > superdiagonal, with the array TAU, represent the unitary */
/* > matrix Q as a product of elementary reflectors;
if UPLO */
/* > = 'L', the diagonal and first subdiagonal of A are over- */
/* > written by the corresponding elements of the tridiagonal */
/* > matrix T, and the elements below the first subdiagonal, with */
/* > the array TAU, represent the unitary matrix Q as a product */
/* > of elementary reflectors. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The diagonal elements of the tridiagonal matrix T: */
/* > D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* > E is REAL array, dimension (N-1) */
/* > The off-diagonal elements of the tridiagonal matrix T: */
/* > E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX array, dimension (N-1) */
/* > The scalar factors of the elementary reflectors (see Further */
/* > Details). */
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
/* > \ingroup complexOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > If UPLO = 'U', the matrix Q is represented as a product of elementary */
/* > reflectors */
/* > */
/* > Q = H(n-1) . . . H(2) H(1). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar, and v is a complex vector with */
/* > v(i+1:n) = 0 and v(i) = 1;
v(1:i-1) is stored on exit in AP, */
/* > overwriting A(1:i-1,i+1), and tau is stored in TAU(i). */
/* > */
/* > If UPLO = 'L', the matrix Q is represented as a product of elementary */
/* > reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(n-1). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar, and v is a complex vector with */
/* > v(1:i) = 0 and v(i+1) = 1;
v(i+2:n) is stored on exit in AP, */
/* > overwriting A(i+2:n,i), and tau is stored in TAU(i). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void chptrd_(char *uplo, integer *n, complex *ap, real *d__, real *e, complex *tau, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("chptrd inputs: uplo %c, n %" FLA_IS "", *uplo, *n);
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;
    complex q__1, q__2, q__3, q__4;
    /* Local variables */
    integer i__, i1, ii, i1i1;
    complex taui;
    extern /* Subroutine */
        void
        chpr2_(char *, integer *, complex *, complex *, integer *, complex *, integer *, complex *);
    complex alpha;
    extern /* Complex */
        VOID
        cdotc_f2c_(complex *, integer *, complex *, integer *, complex *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        chpmv_(char *, integer *, complex *, complex *, complex *, integer *, complex *, complex *,
               integer *),
        caxpy_(integer *, complex *, complex *, integer *, complex *, integer *);
    logical upper;
    extern /* Subroutine */
        void
        clarfg_(integer *, complex *, complex *, integer *, complex *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
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
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    --tau;
    --e;
    --d__;
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
        xerbla_("CHPTRD", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n <= 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(upper)
    {
        /* Reduce the upper triangle of A. */
        /* I1 is the index in AP of A(1,I+1). */
        i1 = *n * (*n - 1) / 2 + 1;
        i__1 = i1 + *n - 1;
        i__2 = i1 + *n - 1;
        r__1 = ap[i__2].r;
        ap[i__1].r = r__1;
        ap[i__1].i = 0.f; // , expr subst
        for(i__ = *n - 1; i__ >= 1; --i__)
        {
            /* Generate elementary reflector H(i) = I - tau * v * v**H */
            /* to annihilate A(1:i-1,i+1) */
            i__1 = i1 + i__ - 1;
            alpha.r = ap[i__1].r;
            alpha.i = ap[i__1].i; // , expr subst
            clarfg_(&i__, &alpha, &ap[i1], &c__1, &taui);
            e[i__] = alpha.r;
            if(taui.r != 0.f || taui.i != 0.f)
            {
                /* Apply H(i) from both sides to A(1:i,1:i) */
                i__1 = i1 + i__ - 1;
                ap[i__1].r = 1.f;
                ap[i__1].i = 0.f; // , expr subst
                /* Compute y := tau * A * v storing y in TAU(1:i) */
                chpmv_(uplo, &i__, &taui, &ap[1], &ap[i1], &c__1, &c_b2, &tau[1], &c__1);
                /* Compute w := y - 1/2 * tau * (y**H *v) * v */
                q__3.r = -.5f;
                q__3.i = -0.f; // , expr subst
                q__2.r = q__3.r * taui.r - q__3.i * taui.i;
                q__2.i = q__3.r * taui.i + q__3.i * taui.r; // , expr subst
                cdotc_f2c_(&q__4, &i__, &tau[1], &c__1, &ap[i1], &c__1);
                q__1.r = q__2.r * q__4.r - q__2.i * q__4.i;
                q__1.i = q__2.r * q__4.i + q__2.i * q__4.r; // , expr subst
                alpha.r = q__1.r;
                alpha.i = q__1.i; // , expr subst
                caxpy_(&i__, &alpha, &ap[i1], &c__1, &tau[1], &c__1);
                /* Apply the transformation as a rank-2 update: */
                /* A := A - v * w**H - w * v**H */
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                chpr2_(uplo, &i__, &q__1, &ap[i1], &c__1, &tau[1], &c__1, &ap[1]);
            }
            i__1 = i1 + i__ - 1;
            i__2 = i__;
            ap[i__1].r = e[i__2];
            ap[i__1].i = 0.f; // , expr subst
            i__1 = i1 + i__;
            d__[i__ + 1] = ap[i__1].r;
            i__1 = i__;
            tau[i__1].r = taui.r;
            tau[i__1].i = taui.i; // , expr subst
            i1 -= i__;
            /* L10: */
        }
        d__[1] = ap[1].r;
    }
    else
    {
        /* Reduce the lower triangle of A. II is the index in AP of */
        /* A(i,i) and I1I1 is the index of A(i+1,i+1). */
        ii = 1;
        r__1 = ap[1].r;
        ap[1].r = r__1;
        ap[1].i = 0.f; // , expr subst
        i__1 = *n - 1;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i1i1 = ii + *n - i__ + 1;
            /* Generate elementary reflector H(i) = I - tau * v * v**H */
            /* to annihilate A(i+2:n,i) */
            i__2 = ii + 1;
            alpha.r = ap[i__2].r;
            alpha.i = ap[i__2].i; // , expr subst
            i__2 = *n - i__;
            clarfg_(&i__2, &alpha, &ap[ii + 2], &c__1, &taui);
            e[i__] = alpha.r;
            if(taui.r != 0.f || taui.i != 0.f)
            {
                /* Apply H(i) from both sides to A(i+1:n,i+1:n) */
                i__2 = ii + 1;
                ap[i__2].r = 1.f;
                ap[i__2].i = 0.f; // , expr subst
                /* Compute y := tau * A * v storing y in TAU(i:n-1) */
                i__2 = *n - i__;
                chpmv_(uplo, &i__2, &taui, &ap[i1i1], &ap[ii + 1], &c__1, &c_b2, &tau[i__], &c__1);
                /* Compute w := y - 1/2 * tau * (y**H *v) * v */
                q__3.r = -.5f;
                q__3.i = -0.f; // , expr subst
                q__2.r = q__3.r * taui.r - q__3.i * taui.i;
                q__2.i = q__3.r * taui.i + q__3.i * taui.r; // , expr subst
                i__2 = *n - i__;
                cdotc_f2c_(&q__4, &i__2, &tau[i__], &c__1, &ap[ii + 1], &c__1);
                q__1.r = q__2.r * q__4.r - q__2.i * q__4.i;
                q__1.i = q__2.r * q__4.i + q__2.i * q__4.r; // , expr subst
                alpha.r = q__1.r;
                alpha.i = q__1.i; // , expr subst
                i__2 = *n - i__;
                caxpy_(&i__2, &alpha, &ap[ii + 1], &c__1, &tau[i__], &c__1);
                /* Apply the transformation as a rank-2 update: */
                /* A := A - v * w**H - w * v**H */
                i__2 = *n - i__;
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                chpr2_(uplo, &i__2, &q__1, &ap[ii + 1], &c__1, &tau[i__], &c__1, &ap[i1i1]);
            }
            i__2 = ii + 1;
            i__3 = i__;
            ap[i__2].r = e[i__3];
            ap[i__2].i = 0.f; // , expr subst
            i__2 = ii;
            d__[i__] = ap[i__2].r;
            i__2 = i__;
            tau[i__2].r = taui.r;
            tau[i__2].i = taui.i; // , expr subst
            ii = i1i1;
            /* L20: */
        }
        i__1 = ii;
        d__[*n] = ap[i__1].r;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of CHPTRD */
}
/* chptrd_ */
