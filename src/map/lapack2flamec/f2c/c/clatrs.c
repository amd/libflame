/* clatrs.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */

/*
 *     Modifications Copyright (c) 2024 Advanced Micro Devices, Inc.  All rights reserved.
 */

#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static real c_b40 = .5f;
/* > \brief \b CLATRS solves a triangular system of equations with the scale factor set to prevent
 * overflow. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLATRS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatrs.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatrs.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatrs.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE, */
/* CNORM, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIAG, NORMIN, TRANS, UPLO */
/* INTEGER INFO, LDA, N */
/* REAL SCALE */
/* .. */
/* .. Array Arguments .. */
/* REAL CNORM( * ) */
/* COMPLEX A( LDA, * ), X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLATRS solves one of the triangular systems */
/* > */
/* > A * x = s*b, A**T * x = s*b, or A**H * x = s*b, */
/* > */
/* > with scaling to prevent overflow. Here A is an upper or lower */
/* > triangular matrix, A**T denotes the transpose of A, A**H denotes the */
/* > conjugate transpose of A, x and b are n-element vectors, and s is a */
/* > scaling factor, usually less than or equal to 1, chosen so that the */
/* > components of x will be less than the overflow threshold. If the */
/* > unscaled problem will not cause overflow, the Level 2 BLAS routine */
/* > CTRSV is called. If the matrix A is singular (A(j,j) = 0 for some j), */
/* > then s is set to 0 and a non-trivial solution to A*x = 0 is returned. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the matrix A is upper or lower triangular. */
/* > = 'U': Upper triangular */
/* > = 'L': Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the operation applied to A. */
/* > = 'N': Solve A * x = s*b (No transpose) */
/* > = 'T': Solve A**T * x = s*b (Transpose) */
/* > = 'C': Solve A**H * x = s*b (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* > DIAG is CHARACTER*1 */
/* > Specifies whether or not the matrix A is unit triangular. */
/* > = 'N': Non-unit triangular */
/* > = 'U': Unit triangular */
/* > \endverbatim */
/* > */
/* > \param[in] NORMIN */
/* > \verbatim */
/* > NORMIN is CHARACTER*1 */
/* > Specifies whether CNORM has been set or not. */
/* > = 'Y': CNORM contains the column norms on entry */
/* > = 'N': CNORM is not set on entry. On exit, the norms will */
/* > be computed and stored in CNORM. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > The triangular matrix A. If UPLO = 'U', the leading n by n */
/* > upper triangular part of the array A contains the upper */
/* > triangular matrix, and the strictly lower triangular part of */
/* > A is not referenced. If UPLO = 'L', the leading n by n lower */
/* > triangular part of the array A contains the lower triangular */
/* > matrix, and the strictly upper triangular part of A is not */
/* > referenced. If DIAG = 'U', the diagonal elements of A are */
/* > also not referenced and are assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max (1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension (N) */
/* > On entry, the right hand side b of the triangular system. */
/* > On exit, X is overwritten by the solution vector x. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* > SCALE is REAL */
/* > The scaling factor s for the triangular system */
/* > A * x = s*b, A**T * x = s*b, or A**H * x = s*b. */
/* > If SCALE = 0, the matrix A is singular or badly scaled, and */
/* > the vector x is an exact or approximate solution to A*x = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CNORM */
/* > \verbatim */
/* > CNORM is REAL array, dimension (N) */
/* > */
/* > If NORMIN = 'Y', CNORM is an input argument and CNORM(j) */
/* > contains the norm of the off-diagonal part of the j-th column */
/* > of A. If TRANS = 'N', CNORM(j) must be greater than or equal */
/* > to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j) */
/* > must be greater than or equal to the 1-norm. */
/* > */
/* > If NORMIN = 'N', CNORM is an output argument and CNORM(j) */
/* > returns the 1-norm of the offdiagonal part of the j-th column */
/* > of A. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -k, the k-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup complexOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > A rough bound on x is computed;
if that is less than overflow, CTRSV */
/* > is called, otherwise, specific code is used which checks for possible */
/* > overflow or divide-by-zero at every operation. */
/* > */
/* > A columnwise scheme is used for solving A*x = b. The basic algorithm */
/* > if A is lower triangular is */
/* > */
/* > x[1:n] := b[1:n] */
/* > for j = 1, ..., n */
/* > x(j) := x(j) / A(j,j) */
/* > x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j] */
/* > end */
/* > */
/* > Define bounds on the components of x after j iterations of the loop: */
/* > M(j) = bound on x[1:j] */
/* > G(j) = bound on x[j+1:n] */
/* > Initially, let M(0) = 0 and G(0) = max{
x(i), i=1,...,n}
. */
/* > */
/* > Then for iteration j+1 we have */
/* > M(j+1) <= G(j) / | A(j+1,j+1) | */
/* > G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] | */
/* > <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | ) */
/* > */
/* > where CNORM(j+1) is greater than or equal to the infinity-norm of */
/* > column j+1 of A, not counting the diagonal. Hence */
/* > */
/* > G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | ) */
/* > 1<=i<=j */
/* > and */
/* > */
/* > |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| ) */
/* > 1<=i< j */
/* > */
/* > Since |x(j)| <= M(j), we use the Level 2 BLAS routine CTRSV if the */
/* > reciprocal of the largest M(j), j=1,..,n, is larger than */
/* > fla_max(underflow, 1/overflow). */
/* > */
/* > The bound on x(j) is also used to determine when a step in the */
/* > columnwise method can be performed without fear of overflow. If */
/* > the computed bound is greater than a large constant, x is scaled to */
/* > prevent overflow, but if the bound overflows, x is set to 0, x(j) to */
/* > 1, and scale to 0, and a non-trivial solution to A*x = 0 is found. */
/* > */
/* > Similarly, a row-wise scheme is used to solve A**T *x = b or */
/* > A**H *x = b. The basic algorithm for A upper triangular is */
/* > */
/* > for j = 1, ..., n */
/* > x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j) */
/* > end */
/* > */
/* > We simultaneously compute two bounds */
/* > G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j */
/* > M(j) = bound on x(i), 1<=i<=j */
/* > */
/* > The initial values are G(0) = 0, M(0) = max{
b(i), i=1,..,n}
, and we */
/* > add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1. */
/* > Then the bound on x(j) is */
/* > */
/* > M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) | */
/* > */
/* > <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| ) */
/* > 1<=i<=j */
/* > */
/* > and we can safely call CTRSV if 1/M(n) and 1/G(n) are both greater */
/* > than fla_max(underflow, 1/overflow). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void clatrs_(char *uplo, char *trans, char *diag, char *normin, integer *n, complex *a,
             integer *lda, complex *x, real *scale, real *cnorm, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clatrs inputs: uplo %c, trans %c, diag %c, normin %c, n %lld, lda %lld",
             *uplo, *trans, *diag, *normin, *n, *lda);
#else
    snprintf(buffer, 256, "clatrs inputs: uplo %c, trans %c, diag %c, normin %c, n %d, lda %d",
             *uplo, *trans, *diag, *normin, *n, *lda);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3, q__4;
    /* Builtin functions */
    double r_imag(complex *);
    void r_cnjg(complex *, complex *);
    /* Local variables */
    integer i__, j;
    real xj, rec, tjj;
    integer jinc;
    real xbnd;
    integer imax;
    real tmax;
    complex tjjs;
    real xmax, grow;
    extern /* Complex */
        VOID
        cdotc_f2c_(complex *, integer *, complex *, integer *, complex *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *);
    real tscal;
    complex uscal;
    integer jlast;
    extern /* Complex */
        VOID
        cdotu_f2c_(complex *, integer *, complex *, integer *, complex *, integer *);
    complex csumj;
    extern /* Subroutine */
        void
        caxpy_(integer *, complex *, complex *, integer *, complex *, integer *);
    logical upper;
    extern /* Subroutine */
        void
        ctrsv_(char *, char *, char *, integer *, complex *, integer *, complex *, integer *);
    extern integer icamax_(integer *, complex *, integer *);
    extern /* Complex */
        void
        cladiv_f2c_(complex *, complex *, complex *);
    extern real slamch_(char *);
    extern /* Subroutine */
        void
        csscal_(integer *, real *, complex *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real bignum;
    extern integer isamax_(integer *, real *, integer *);
    extern real scasum_(integer *, complex *, integer *);
    logical notran;
    integer jfirst;
    real smlnum;
    logical nounit;
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
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --cnorm;
    /* Function Body */
    *info = 0;
    // initializing as {1, 0} because it is
    // used as divisor
    tjjs = (complex){.r = 1.f, .i = 0.f};
    upper = lsame_(uplo, "U", 1, 1);
    notran = lsame_(trans, "N", 1, 1);
    nounit = lsame_(diag, "N", 1, 1);
    /* Test the input parameters. */
    if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(!notran && !lsame_(trans, "T", 1, 1) && !lsame_(trans, "C", 1, 1))
    {
        *info = -2;
    }
    else if(!nounit && !lsame_(diag, "U", 1, 1))
    {
        *info = -3;
    }
    else if(!lsame_(normin, "Y", 1, 1) && !lsame_(normin, "N", 1, 1))
    {
        *info = -4;
    }
    else if(*n < 0)
    {
        *info = -5;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -7;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CLATRS", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    *scale = 1.f;
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Determine machine dependent parameters to control overflow. */
    smlnum = slamch_("Safe minimum") / slamch_("Precision");
    bignum = 1.f / smlnum;
    if(lsame_(normin, "N", 1, 1))
    {
        /* Compute the 1-norm of each column, not including the diagonal. */
        if(upper)
        {
            /* A is upper triangular. */
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = j - 1;
                cnorm[j] = scasum_(&i__2, &a[j * a_dim1 + 1], &c__1);
                /* L10: */
            }
        }
        else
        {
            /* A is lower triangular. */
            i__1 = *n - 1;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = *n - j;
                cnorm[j] = scasum_(&i__2, &a[j + 1 + j * a_dim1], &c__1);
                /* L20: */
            }
            cnorm[*n] = 0.f;
        }
    }
    /* Scale the column norms by TSCAL if the maximum element in CNORM is */
    /* greater than BIGNUM/2. */
    imax = isamax_(n, &cnorm[1], &c__1);
    tmax = cnorm[imax];
    if(tmax <= bignum * .5f)
    {
        tscal = 1.f;
    }
    else
    {
        /* Avoid NaN generation if entries in CNORM exceed the */
        /* overflow threshold */
        if(tmax <= slamch_("Overflow"))
        {
            /* Case 1: All entries in CNORM are valid floating-point numbers */
            tscal = .5f / (smlnum * tmax);
            sscal_(n, &tscal, &cnorm[1], &c__1);
        }
        else
        {
            /* Case 2: At least one column norm of A cannot be */
            /* represented as a floating-point number. Find the */
            /* maximum offdiagonal absolute value */
            /* fla_max( |Re(A(I,J))|, |Im(A(I,J)| ). If this entry is */
            /* not +/- Infinity, use this value as TSCAL. */
            tmax = 0.f;
            if(upper)
            {
                /* A is upper triangular. */
                i__1 = *n;
                for(j = 2; j <= i__1; ++j)
                {
                    i__2 = j - 1;
                    for(i__ = 1; i__ <= i__2; ++i__)
                    {
                        /* Computing MAX */
                        i__3 = i__ + j * a_dim1;
                        r__3 = tmax, r__4 = (r__1 = a[i__3].r, f2c_abs(r__1));
                        r__3 = fla_max(r__3, r__4);
                        r__4 = (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2)); // ; expr subst
                        tmax = fla_max(r__3, r__4);
                    }
                }
            }
            else
            {
                /* A is lower triangular. */
                i__1 = *n - 1;
                for(j = 1; j <= i__1; ++j)
                {
                    i__2 = *n;
                    for(i__ = j + 1; i__ <= i__2; ++i__)
                    {
                        /* Computing MAX */
                        i__3 = i__ + j * a_dim1;
                        r__3 = tmax, r__4 = (r__1 = a[i__3].r, f2c_abs(r__1));
                        r__3 = fla_max(r__3, r__4);
                        r__4 = (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2)); // ; expr subst
                        tmax = fla_max(r__3, r__4);
                    }
                }
            }
            if(tmax <= slamch_("Overflow"))
            {
                tscal = 1.f / (smlnum * tmax);
                i__1 = *n;
                for(j = 1; j <= i__1; ++j)
                {
                    if(cnorm[j] <= slamch_("Overflow"))
                    {
                        cnorm[j] *= tscal;
                    }
                    else
                    {
                        /* Recompute the 1-norm of each column without */
                        /* introducing Infinity in the summation. */
                        tscal *= 2.f;
                        cnorm[j] = 0.f;
                        if(upper)
                        {
                            i__2 = j - 1;
                            for(i__ = 1; i__ <= i__2; ++i__)
                            {
                                i__3 = i__ + j * a_dim1;
                                cnorm[j] += tscal
                                            * ((r__1 = a[i__3].r / 2.f, f2c_abs(r__1))
                                               + (r__2 = r_imag(&a[i__ + j * a_dim1]) / 2.f,
                                                  f2c_abs(r__2)));
                            }
                        }
                        else
                        {
                            i__2 = *n;
                            for(i__ = j + 1; i__ <= i__2; ++i__)
                            {
                                i__3 = i__ + j * a_dim1;
                                cnorm[j] += tscal
                                            * ((r__1 = a[i__3].r / 2.f, f2c_abs(r__1))
                                               + (r__2 = r_imag(&a[i__ + j * a_dim1]) / 2.f,
                                                  f2c_abs(r__2)));
                            }
                        }
                        tscal *= .5f;
                    }
                }
            }
            else
            {
                /* At least one entry of A is not a valid floating-point */
                /* entry. Rely on TRSV to propagate Inf and NaN. */
                ctrsv_(uplo, trans, diag, n, &a[a_offset], lda, &x[1], &c__1);
                AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
                return;
            }
        }
    }
    /* Compute a bound on the computed solution vector to see if the */
    /* Level 2 BLAS routine CTRSV can be used. */
    xmax = 0.f;
    i__1 = *n;
    for(j = 1; j <= i__1; ++j)
    {
        /* Computing MAX */
        i__2 = j;
        r__3 = xmax;
        r__4 = (r__1 = x[i__2].r / 2.f, f2c_abs(r__1))
               + (r__2 = x[j].i / 2.f, f2c_abs(r__2)); // , expr subst
        xmax = fla_max(r__3, r__4);
        /* L30: */
    }
    xbnd = xmax;
    if(notran)
    {
        /* Compute the growth in A * x = b. */
        if(upper)
        {
            jfirst = *n;
            jlast = 1;
            jinc = -1;
        }
        else
        {
            jfirst = 1;
            jlast = *n;
            jinc = 1;
        }
        if(tscal != 1.f)
        {
            grow = 0.f;
            goto L60;
        }
        if(nounit)
        {
            /* A is non-unit triangular. */
            /* Compute GROW = 1/G(j) and XBND = 1/M(j). */
            /* Initially, G(0) = max{
            x(i), i=1,...,n}
            . */
            grow = .5f / fla_max(xbnd, smlnum);
            xbnd = grow;
            i__1 = jlast;
            i__2 = jinc;
            for(j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
            {
                /* Exit the loop if the growth factor is too small. */
                if(grow <= smlnum)
                {
                    goto L60;
                }
                i__3 = j + j * a_dim1;
                tjjs.r = a[i__3].r;
                tjjs.i = a[i__3].i; // , expr subst
                tjj = (r__1 = tjjs.r, f2c_abs(r__1)) + (r__2 = tjjs.i, f2c_abs(r__2));
                if(tjj >= smlnum)
                {
                    /* M(j) = G(j-1) / f2c_abs(A(j,j)) */
                    /* Computing MIN */
                    r__1 = xbnd;
                    r__2 = fla_min(1.f, tjj) * grow; // , expr subst
                    xbnd = fla_min(r__1, r__2);
                }
                else
                {
                    /* M(j) could overflow, set XBND to 0. */
                    xbnd = 0.f;
                }
                if(tjj + cnorm[j] >= smlnum)
                {
                    /* G(j) = G(j-1)*( 1 + CNORM(j) / f2c_abs(A(j,j)) ) */
                    grow *= tjj / (tjj + cnorm[j]);
                }
                else
                {
                    /* G(j) could overflow, set GROW to 0. */
                    grow = 0.f;
                }
                /* L40: */
            }
            grow = xbnd;
        }
        else
        {
            /* A is unit triangular. */
            /* Compute GROW = 1/G(j), where G(0) = max{
            x(i), i=1,...,n}
            . */
            /* Computing MIN */
            r__1 = 1.f;
            r__2 = .5f / fla_max(xbnd, smlnum); // , expr subst
            grow = fla_min(r__1, r__2);
            i__2 = jlast;
            i__1 = jinc;
            for(j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1)
            {
                /* Exit the loop if the growth factor is too small. */
                if(grow <= smlnum)
                {
                    goto L60;
                }
                /* G(j) = G(j-1)*( 1 + CNORM(j) ) */
                grow *= 1.f / (cnorm[j] + 1.f);
                /* L50: */
            }
        }
    L60:;
    }
    else
    {
        /* Compute the growth in A**T * x = b or A**H * x = b. */
        if(upper)
        {
            jfirst = 1;
            jlast = *n;
            jinc = 1;
        }
        else
        {
            jfirst = *n;
            jlast = 1;
            jinc = -1;
        }
        if(tscal != 1.f)
        {
            grow = 0.f;
            goto L90;
        }
        if(nounit)
        {
            /* A is non-unit triangular. */
            /* Compute GROW = 1/G(j) and XBND = 1/M(j). */
            /* Initially, M(0) = max{
            x(i), i=1,...,n}
            . */
            grow = .5f / fla_max(xbnd, smlnum);
            xbnd = grow;
            i__1 = jlast;
            i__2 = jinc;
            for(j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
            {
                /* Exit the loop if the growth factor is too small. */
                if(grow <= smlnum)
                {
                    goto L90;
                }
                /* G(j) = fla_max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */
                xj = cnorm[j] + 1.f;
                /* Computing MIN */
                r__1 = grow;
                r__2 = xbnd / xj; // , expr subst
                grow = fla_min(r__1, r__2);
                i__3 = j + j * a_dim1;
                tjjs.r = a[i__3].r;
                tjjs.i = a[i__3].i; // , expr subst
                tjj = (r__1 = tjjs.r, f2c_abs(r__1)) + (r__2 = tjjs.i, f2c_abs(r__2));
                if(tjj >= smlnum)
                {
                    /* M(j) = M(j-1)*( 1 + CNORM(j) ) / f2c_abs(A(j,j)) */
                    if(xj > tjj)
                    {
                        xbnd *= tjj / xj;
                    }
                }
                else
                {
                    /* M(j) could overflow, set XBND to 0. */
                    xbnd = 0.f;
                }
                /* L70: */
            }
            grow = fla_min(grow, xbnd);
        }
        else
        {
            /* A is unit triangular. */
            /* Compute GROW = 1/G(j), where G(0) = max{
            x(i), i=1,...,n}
            . */
            /* Computing MIN */
            r__1 = 1.f;
            r__2 = .5f / fla_max(xbnd, smlnum); // , expr subst
            grow = fla_min(r__1, r__2);
            i__2 = jlast;
            i__1 = jinc;
            for(j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1)
            {
                /* Exit the loop if the growth factor is too small. */
                if(grow <= smlnum)
                {
                    goto L90;
                }
                /* G(j) = ( 1 + CNORM(j) )*G(j-1) */
                xj = cnorm[j] + 1.f;
                grow /= xj;
                /* L80: */
            }
        }
    L90:;
    }
    if(grow * tscal > smlnum)
    {
        /* Use the Level 2 BLAS solve if the reciprocal of the bound on */
        /* elements of X is not too small. */
        ctrsv_(uplo, trans, diag, n, &a[a_offset], lda, &x[1], &c__1);
    }
    else
    {
        /* Use a Level 1 BLAS solve, scaling intermediate results. */
        if(xmax > bignum * .5f)
        {
            /* Scale X so that its components are less than or equal to */
            /* BIGNUM in absolute value. */
            *scale = bignum * .5f / xmax;
            csscal_(n, scale, &x[1], &c__1);
            xmax = bignum;
        }
        else
        {
            xmax *= 2.f;
        }
        if(notran)
        {
            /* Solve A * x = b */
            i__1 = jlast;
            i__2 = jinc;
            for(j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
            {
                /* Compute x(j) = b(j) / A(j,j), scaling x if necessary. */
                i__3 = j;
                xj = (r__1 = x[i__3].r, f2c_abs(r__1)) + (r__2 = x[j].i, f2c_abs(r__2));
                if(nounit)
                {
                    i__3 = j + j * a_dim1;
                    q__1.r = tscal * a[i__3].r;
                    q__1.i = tscal * a[i__3].i; // , expr subst
                    tjjs.r = q__1.r;
                    tjjs.i = q__1.i; // , expr subst
                }
                else
                {
                    tjjs.r = tscal;
                    tjjs.i = 0.f; // , expr subst
                    if(tscal == 1.f)
                    {
                        goto L105;
                    }
                }
                tjj = (r__1 = tjjs.r, f2c_abs(r__1)) + (r__2 = tjjs.i, f2c_abs(r__2));
                if(tjj > smlnum)
                {
                    /* f2c_abs(A(j,j)) > SMLNUM: */
                    if(tjj < 1.f)
                    {
                        if(xj > tjj * bignum)
                        {
                            /* Scale x by 1/b(j). */
                            rec = 1.f / xj;
                            csscal_(n, &rec, &x[1], &c__1);
                            *scale *= rec;
                            xmax *= rec;
                        }
                    }
                    i__3 = j;
                    cladiv_f2c_(&q__1, &x[j], &tjjs);
                    x[i__3].r = q__1.r;
                    x[i__3].i = q__1.i; // , expr subst
                    i__3 = j;
                    xj = (r__1 = x[i__3].r, f2c_abs(r__1)) + (r__2 = x[j].i, f2c_abs(r__2));
                }
                else if(tjj > 0.f)
                {
                    /* 0 < f2c_abs(A(j,j)) <= SMLNUM: */
                    if(xj > tjj * bignum)
                    {
                        /* Scale x by (1/f2c_abs(x(j)))*f2c_abs(A(j,j))*BIGNUM */
                        /* to avoid overflow when dividing by A(j,j). */
                        rec = tjj * bignum / xj;
                        if(cnorm[j] > 1.f)
                        {
                            /* Scale by 1/CNORM(j) to avoid overflow when */
                            /* multiplying x(j) times column j. */
                            rec /= cnorm[j];
                        }
                        csscal_(n, &rec, &x[1], &c__1);
                        *scale *= rec;
                        xmax *= rec;
                    }
                    i__3 = j;
                    cladiv_f2c_(&q__1, &x[j], &tjjs);
                    x[i__3].r = q__1.r;
                    x[i__3].i = q__1.i; // , expr subst
                    i__3 = j;
                    xj = (r__1 = x[i__3].r, f2c_abs(r__1)) + (r__2 = x[j].i, f2c_abs(r__2));
                }
                else
                {
                    /* A(j,j) = 0: Set x(1:n) = 0, x(j) = 1, and */
                    /* scale = 0, and compute a solution to A*x = 0. */
                    i__3 = *n;
                    for(i__ = 1; i__ <= i__3; ++i__)
                    {
                        i__4 = i__;
                        x[i__4].r = 0.f;
                        x[i__4].i = 0.f; // , expr subst
                        /* L100: */
                    }
                    i__3 = j;
                    x[i__3].r = 1.f;
                    x[i__3].i = 0.f; // , expr subst
                    xj = 1.f;
                    *scale = 0.f;
                    xmax = 0.f;
                }
            L105: /* Scale x if necessary to avoid overflow when adding a */
                /* multiple of column j of A. */
                if(xj > 1.f)
                {
                    rec = 1.f / xj;
                    if(cnorm[j] > (bignum - xmax) * rec)
                    {
                        /* Scale x by 1/(2*f2c_abs(x(j))). */
                        rec *= .5f;
                        csscal_(n, &rec, &x[1], &c__1);
                        *scale *= rec;
                    }
                }
                else if(xj * cnorm[j] > bignum - xmax)
                {
                    /* Scale x by 1/2. */
                    csscal_(n, &c_b40, &x[1], &c__1);
                    *scale *= .5f;
                }
                if(upper)
                {
                    if(j > 1)
                    {
                        /* Compute the update */
                        /* x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j) */
                        i__3 = j - 1;
                        i__4 = j;
                        q__2.r = -x[i__4].r;
                        q__2.i = -x[i__4].i; // , expr subst
                        q__1.r = tscal * q__2.r;
                        q__1.i = tscal * q__2.i; // , expr subst
                        caxpy_(&i__3, &q__1, &a[j * a_dim1 + 1], &c__1, &x[1], &c__1);
                        i__3 = j - 1;
                        i__ = icamax_(&i__3, &x[1], &c__1);
                        i__3 = i__;
                        xmax = (r__1 = x[i__3].r, f2c_abs(r__1)) + (r__2 = x[i__].i, f2c_abs(r__2));
                    }
                }
                else
                {
                    if(j < *n)
                    {
                        /* Compute the update */
                        /* x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j) */
                        i__3 = *n - j;
                        i__4 = j;
                        q__2.r = -x[i__4].r;
                        q__2.i = -x[i__4].i; // , expr subst
                        q__1.r = tscal * q__2.r;
                        q__1.i = tscal * q__2.i; // , expr subst
                        caxpy_(&i__3, &q__1, &a[j + 1 + j * a_dim1], &c__1, &x[j + 1], &c__1);
                        i__3 = *n - j;
                        i__ = j + icamax_(&i__3, &x[j + 1], &c__1);
                        i__3 = i__;
                        xmax = (r__1 = x[i__3].r, f2c_abs(r__1)) + (r__2 = x[i__].i, f2c_abs(r__2));
                    }
                }
                /* L110: */
            }
        }
        else if(lsame_(trans, "T", 1, 1))
        {
            /* Solve A**T * x = b */
            i__2 = jlast;
            i__1 = jinc;
            for(j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1)
            {
                /* Compute x(j) = b(j) - sum A(k,j)*x(k). */
                /* k<>j */
                i__3 = j;
                xj = (r__1 = x[i__3].r, f2c_abs(r__1)) + (r__2 = x[j].i, f2c_abs(r__2));
                uscal.r = tscal;
                uscal.i = 0.f; // , expr subst
                rec = 1.f / fla_max(xmax, 1.f);
                if(cnorm[j] > (bignum - xj) * rec)
                {
                    /* If x(j) could overflow, scale x by 1/(2*XMAX). */
                    rec *= .5f;
                    if(nounit)
                    {
                        i__3 = j + j * a_dim1;
                        q__1.r = tscal * a[i__3].r;
                        q__1.i = tscal * a[i__3].i; // , expr subst
                        tjjs.r = q__1.r;
                        tjjs.i = q__1.i; // , expr subst
                    }
                    else
                    {
                        tjjs.r = tscal;
                        tjjs.i = 0.f; // , expr subst
                    }
                    tjj = (r__1 = tjjs.r, f2c_abs(r__1)) + (r__2 = tjjs.i, f2c_abs(r__2));
                    if(tjj > 1.f)
                    {
                        /* Divide by A(j,j) when scaling x if A(j,j) > 1. */
                        /* Computing MIN */
                        r__1 = 1.f;
                        r__2 = rec * tjj; // , expr subst
                        rec = fla_min(r__1, r__2);
                        cladiv_f2c_(&q__1, &uscal, &tjjs);
                        uscal.r = q__1.r;
                        uscal.i = q__1.i; // , expr subst
                    }
                    if(rec < 1.f)
                    {
                        csscal_(n, &rec, &x[1], &c__1);
                        *scale *= rec;
                        xmax *= rec;
                    }
                }
                csumj.r = 0.f;
                csumj.i = 0.f; // , expr subst
                if(uscal.r == 1.f && uscal.i == 0.f)
                {
                    /* If the scaling needed for A in the dot product is 1, */
                    /* call CDOTU to perform the dot product. */
                    if(upper)
                    {
                        i__3 = j - 1;
                        cdotu_f2c_(&q__1, &i__3, &a[j * a_dim1 + 1], &c__1, &x[1], &c__1);
                        csumj.r = q__1.r;
                        csumj.i = q__1.i; // , expr subst
                    }
                    else if(j < *n)
                    {
                        i__3 = *n - j;
                        cdotu_f2c_(&q__1, &i__3, &a[j + 1 + j * a_dim1], &c__1, &x[j + 1], &c__1);
                        csumj.r = q__1.r;
                        csumj.i = q__1.i; // , expr subst
                    }
                }
                else
                {
                    /* Otherwise, use in-line code for the dot product. */
                    if(upper)
                    {
                        i__3 = j - 1;
                        for(i__ = 1; i__ <= i__3; ++i__)
                        {
                            i__4 = i__ + j * a_dim1;
                            q__3.r = a[i__4].r * uscal.r - a[i__4].i * uscal.i;
                            q__3.i = a[i__4].r * uscal.i + a[i__4].i * uscal.r; // , expr subst
                            i__5 = i__;
                            q__2.r = q__3.r * x[i__5].r - q__3.i * x[i__5].i;
                            q__2.i = q__3.r * x[i__5].i + q__3.i * x[i__5].r; // , expr subst
                            q__1.r = csumj.r + q__2.r;
                            q__1.i = csumj.i + q__2.i; // , expr subst
                            csumj.r = q__1.r;
                            csumj.i = q__1.i; // , expr subst
                            /* L120: */
                        }
                    }
                    else if(j < *n)
                    {
                        i__3 = *n;
                        for(i__ = j + 1; i__ <= i__3; ++i__)
                        {
                            i__4 = i__ + j * a_dim1;
                            q__3.r = a[i__4].r * uscal.r - a[i__4].i * uscal.i;
                            q__3.i = a[i__4].r * uscal.i + a[i__4].i * uscal.r; // , expr subst
                            i__5 = i__;
                            q__2.r = q__3.r * x[i__5].r - q__3.i * x[i__5].i;
                            q__2.i = q__3.r * x[i__5].i + q__3.i * x[i__5].r; // , expr subst
                            q__1.r = csumj.r + q__2.r;
                            q__1.i = csumj.i + q__2.i; // , expr subst
                            csumj.r = q__1.r;
                            csumj.i = q__1.i; // , expr subst
                            /* L130: */
                        }
                    }
                }
                q__1.r = tscal;
                q__1.i = 0.f; // , expr subst
                if(uscal.r == q__1.r && uscal.i == q__1.i)
                {
                    /* Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j) */
                    /* was not used to scale the dotproduct. */
                    i__3 = j;
                    i__4 = j;
                    q__1.r = x[i__4].r - csumj.r;
                    q__1.i = x[i__4].i - csumj.i; // , expr subst
                    x[i__3].r = q__1.r;
                    x[i__3].i = q__1.i; // , expr subst
                    i__3 = j;
                    xj = (r__1 = x[i__3].r, f2c_abs(r__1)) + (r__2 = x[j].i, f2c_abs(r__2));
                    if(nounit)
                    {
                        i__3 = j + j * a_dim1;
                        q__1.r = tscal * a[i__3].r;
                        q__1.i = tscal * a[i__3].i; // , expr subst
                        tjjs.r = q__1.r;
                        tjjs.i = q__1.i; // , expr subst
                    }
                    else
                    {
                        tjjs.r = tscal;
                        tjjs.i = 0.f; // , expr subst
                        if(tscal == 1.f)
                        {
                            goto L145;
                        }
                    }
                    /* Compute x(j) = x(j) / A(j,j), scaling if necessary. */
                    tjj = (r__1 = tjjs.r, f2c_abs(r__1)) + (r__2 = tjjs.i, f2c_abs(r__2));
                    if(tjj > smlnum)
                    {
                        /* f2c_abs(A(j,j)) > SMLNUM: */
                        if(tjj < 1.f)
                        {
                            if(xj > tjj * bignum)
                            {
                                /* Scale X by 1/f2c_abs(x(j)). */
                                rec = 1.f / xj;
                                csscal_(n, &rec, &x[1], &c__1);
                                *scale *= rec;
                                xmax *= rec;
                            }
                        }
                        i__3 = j;
                        cladiv_f2c_(&q__1, &x[j], &tjjs);
                        x[i__3].r = q__1.r;
                        x[i__3].i = q__1.i; // , expr subst
                    }
                    else if(tjj > 0.f)
                    {
                        /* 0 < f2c_abs(A(j,j)) <= SMLNUM: */
                        if(xj > tjj * bignum)
                        {
                            /* Scale x by (1/f2c_abs(x(j)))*f2c_abs(A(j,j))*BIGNUM. */
                            rec = tjj * bignum / xj;
                            csscal_(n, &rec, &x[1], &c__1);
                            *scale *= rec;
                            xmax *= rec;
                        }
                        i__3 = j;
                        cladiv_f2c_(&q__1, &x[j], &tjjs);
                        x[i__3].r = q__1.r;
                        x[i__3].i = q__1.i; // , expr subst
                    }
                    else
                    {
                        /* A(j,j) = 0: Set x(1:n) = 0, x(j) = 1, and */
                        /* scale = 0 and compute a solution to A**T *x = 0. */
                        i__3 = *n;
                        for(i__ = 1; i__ <= i__3; ++i__)
                        {
                            i__4 = i__;
                            x[i__4].r = 0.f;
                            x[i__4].i = 0.f; // , expr subst
                            /* L140: */
                        }
                        i__3 = j;
                        x[i__3].r = 1.f;
                        x[i__3].i = 0.f; // , expr subst
                        *scale = 0.f;
                        xmax = 0.f;
                    }
                L145:;
                }
                else
                {
                    /* Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot */
                    /* product has already been divided by 1/A(j,j). */
                    i__3 = j;
                    cladiv_f2c_(&q__2, &x[j], &tjjs);
                    q__1.r = q__2.r - csumj.r;
                    q__1.i = q__2.i - csumj.i; // , expr subst
                    x[i__3].r = q__1.r;
                    x[i__3].i = q__1.i; // , expr subst
                }
                /* Computing MAX */
                i__3 = j;
                r__3 = xmax;
                r__4 = (r__1 = x[i__3].r, f2c_abs(r__1))
                       + (r__2 = x[j].i, f2c_abs(r__2)); // , expr subst
                xmax = fla_max(r__3, r__4);
                /* L150: */
            }
        }
        else
        {
            /* Solve A**H * x = b */
            i__1 = jlast;
            i__2 = jinc;
            for(j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
            {
                /* Compute x(j) = b(j) - sum A(k,j)*x(k). */
                /* k<>j */
                i__3 = j;
                xj = (r__1 = x[i__3].r, f2c_abs(r__1)) + (r__2 = x[j].i, f2c_abs(r__2));
                uscal.r = tscal;
                uscal.i = 0.f; // , expr subst
                rec = 1.f / fla_max(xmax, 1.f);
                if(cnorm[j] > (bignum - xj) * rec)
                {
                    /* If x(j) could overflow, scale x by 1/(2*XMAX). */
                    rec *= .5f;
                    if(nounit)
                    {
                        q__2.r = a[j + j * a_dim1].r;
                        q__2.i = -a[j + j * a_dim1].i;
                        q__1.r = tscal * q__2.r;
                        q__1.i = tscal * q__2.i; // , expr subst
                        tjjs.r = q__1.r;
                        tjjs.i = q__1.i; // , expr subst
                    }
                    else
                    {
                        tjjs.r = tscal;
                        tjjs.i = 0.f; // , expr subst
                    }
                    tjj = (r__1 = tjjs.r, f2c_abs(r__1)) + (r__2 = tjjs.i, f2c_abs(r__2));
                    if(tjj > 1.f)
                    {
                        /* Divide by A(j,j) when scaling x if A(j,j) > 1. */
                        /* Computing MIN */
                        r__1 = 1.f;
                        r__2 = rec * tjj; // , expr subst
                        rec = fla_min(r__1, r__2);
                        cladiv_f2c_(&q__1, &uscal, &tjjs);
                        uscal.r = q__1.r;
                        uscal.i = q__1.i; // , expr subst
                    }
                    if(rec < 1.f)
                    {
                        csscal_(n, &rec, &x[1], &c__1);
                        *scale *= rec;
                        xmax *= rec;
                    }
                }
                csumj.r = 0.f;
                csumj.i = 0.f; // , expr subst
                if(uscal.r == 1.f && uscal.i == 0.f)
                {
                    /* If the scaling needed for A in the dot product is 1, */
                    /* call CDOTC to perform the dot product. */
                    if(upper)
                    {
                        i__3 = j - 1;
                        cdotc_f2c_(&q__1, &i__3, &a[j * a_dim1 + 1], &c__1, &x[1], &c__1);
                        csumj.r = q__1.r;
                        csumj.i = q__1.i; // , expr subst
                    }
                    else if(j < *n)
                    {
                        i__3 = *n - j;
                        cdotc_f2c_(&q__1, &i__3, &a[j + 1 + j * a_dim1], &c__1, &x[j + 1], &c__1);
                        csumj.r = q__1.r;
                        csumj.i = q__1.i; // , expr subst
                    }
                }
                else
                {
                    /* Otherwise, use in-line code for the dot product. */
                    if(upper)
                    {
                        i__3 = j - 1;
                        for(i__ = 1; i__ <= i__3; ++i__)
                        {
                            q__4.r = a[i__ + j * a_dim1].r;
                            q__4.i = -a[i__ + j * a_dim1].i;
                            q__3.r = q__4.r * uscal.r - q__4.i * uscal.i;
                            q__3.i = q__4.r * uscal.i + q__4.i * uscal.r; // , expr subst
                            i__4 = i__;
                            q__2.r = q__3.r * x[i__4].r - q__3.i * x[i__4].i;
                            q__2.i = q__3.r * x[i__4].i + q__3.i * x[i__4].r; // , expr subst
                            q__1.r = csumj.r + q__2.r;
                            q__1.i = csumj.i + q__2.i; // , expr subst
                            csumj.r = q__1.r;
                            csumj.i = q__1.i; // , expr subst
                            /* L160: */
                        }
                    }
                    else if(j < *n)
                    {
                        i__3 = *n;
                        for(i__ = j + 1; i__ <= i__3; ++i__)
                        {
                            q__4.r = a[i__ + j * a_dim1].r;
                            q__4.i = -a[i__ + j * a_dim1].i;
                            q__3.r = q__4.r * uscal.r - q__4.i * uscal.i;
                            q__3.i = q__4.r * uscal.i + q__4.i * uscal.r; // , expr subst
                            i__4 = i__;
                            q__2.r = q__3.r * x[i__4].r - q__3.i * x[i__4].i;
                            q__2.i = q__3.r * x[i__4].i + q__3.i * x[i__4].r; // , expr subst
                            q__1.r = csumj.r + q__2.r;
                            q__1.i = csumj.i + q__2.i; // , expr subst
                            csumj.r = q__1.r;
                            csumj.i = q__1.i; // , expr subst
                            /* L170: */
                        }
                    }
                }
                q__1.r = tscal;
                q__1.i = 0.f; // , expr subst
                if(uscal.r == q__1.r && uscal.i == q__1.i)
                {
                    /* Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j) */
                    /* was not used to scale the dotproduct. */
                    i__3 = j;
                    i__4 = j;
                    q__1.r = x[i__4].r - csumj.r;
                    q__1.i = x[i__4].i - csumj.i; // , expr subst
                    x[i__3].r = q__1.r;
                    x[i__3].i = q__1.i; // , expr subst
                    i__3 = j;
                    xj = (r__1 = x[i__3].r, f2c_abs(r__1)) + (r__2 = x[j].i, f2c_abs(r__2));
                    if(nounit)
                    {
                        q__2.r = a[j + j * a_dim1].r;
                        q__2.i = -a[j + j * a_dim1].i;
                        q__1.r = tscal * q__2.r;
                        q__1.i = tscal * q__2.i; // , expr subst
                        tjjs.r = q__1.r;
                        tjjs.i = q__1.i; // , expr subst
                    }
                    else
                    {
                        tjjs.r = tscal;
                        tjjs.i = 0.f; // , expr subst
                        if(tscal == 1.f)
                        {
                            goto L185;
                        }
                    }
                    /* Compute x(j) = x(j) / A(j,j), scaling if necessary. */
                    tjj = (r__1 = tjjs.r, f2c_abs(r__1)) + (r__2 = tjjs.i, f2c_abs(r__2));
                    if(tjj > smlnum)
                    {
                        /* f2c_abs(A(j,j)) > SMLNUM: */
                        if(tjj < 1.f)
                        {
                            if(xj > tjj * bignum)
                            {
                                /* Scale X by 1/f2c_abs(x(j)). */
                                rec = 1.f / xj;
                                csscal_(n, &rec, &x[1], &c__1);
                                *scale *= rec;
                                xmax *= rec;
                            }
                        }
                        i__3 = j;
                        cladiv_f2c_(&q__1, &x[j], &tjjs);
                        x[i__3].r = q__1.r;
                        x[i__3].i = q__1.i; // , expr subst
                    }
                    else if(tjj > 0.f)
                    {
                        /* 0 < f2c_abs(A(j,j)) <= SMLNUM: */
                        if(xj > tjj * bignum)
                        {
                            /* Scale x by (1/f2c_abs(x(j)))*f2c_abs(A(j,j))*BIGNUM. */
                            rec = tjj * bignum / xj;
                            csscal_(n, &rec, &x[1], &c__1);
                            *scale *= rec;
                            xmax *= rec;
                        }
                        i__3 = j;
                        cladiv_f2c_(&q__1, &x[j], &tjjs);
                        x[i__3].r = q__1.r;
                        x[i__3].i = q__1.i; // , expr subst
                    }
                    else
                    {
                        /* A(j,j) = 0: Set x(1:n) = 0, x(j) = 1, and */
                        /* scale = 0 and compute a solution to A**H *x = 0. */
                        i__3 = *n;
                        for(i__ = 1; i__ <= i__3; ++i__)
                        {
                            i__4 = i__;
                            x[i__4].r = 0.f;
                            x[i__4].i = 0.f; // , expr subst
                            /* L180: */
                        }
                        i__3 = j;
                        x[i__3].r = 1.f;
                        x[i__3].i = 0.f; // , expr subst
                        *scale = 0.f;
                        xmax = 0.f;
                    }
                L185:;
                }
                else
                {
                    /* Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot */
                    /* product has already been divided by 1/A(j,j). */
                    i__3 = j;
                    cladiv_f2c_(&q__2, &x[j], &tjjs);
                    q__1.r = q__2.r - csumj.r;
                    q__1.i = q__2.i - csumj.i; // , expr subst
                    x[i__3].r = q__1.r;
                    x[i__3].i = q__1.i; // , expr subst
                }
                /* Computing MAX */
                i__3 = j;
                r__3 = xmax;
                r__4 = (r__1 = x[i__3].r, f2c_abs(r__1))
                       + (r__2 = x[j].i, f2c_abs(r__2)); // , expr subst
                xmax = fla_max(r__3, r__4);
                /* L190: */
            }
        }
        *scale /= tscal;
    }
    /* Scale the column norms by 1/TSCAL for return. */
    if(tscal != 1.f)
    {
        r__1 = 1.f / tscal;
        sscal_(n, &r__1, &cnorm[1], &c__1);
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLATRS */
}
/* clatrs_ */
