/* ./zlatps.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static doublereal c_b36 = .5;
/* > \brief \b ZLATPS solves a triangular system of equations with the matrix held in packed
 * storage. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLATPS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlatps.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlatps.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlatps.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLATPS( UPLO, TRANS, DIAG, NORMIN, N, AP, X, SCALE, */
/* CNORM, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIAG, NORMIN, TRANS, UPLO */
/* INTEGER INFO, N */
/* DOUBLE PRECISION SCALE */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION CNORM( * ) */
/* COMPLEX*16 AP( * ), X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLATPS solves one of the triangular systems */
/* > */
/* > A * x = s*b, A**T * x = s*b, or A**H * x = s*b, */
/* > */
/* > with scaling to prevent overflow, where A is an upper or lower */
/* > triangular matrix stored in packed form. Here A**T denotes the */
/* > transpose of A, A**H denotes the conjugate transpose of A, x and b */
/* > are n-element vectors, and s is a scaling factor, usually less than */
/* > or equal to 1, chosen so that the components of x will be less than */
/* > the overflow threshold. If the unscaled problem will not cause */
/* > overflow, the Level 2 BLAS routine ZTPSV is called. If the matrix A */
/* > is singular (A(j,j) = 0 for some j), then s is set to 0 and a */
/* > non-trivial solution to A*x = 0 is returned. */
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
/* > \param[in] AP */
/* > \verbatim */
/* > AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* > The upper or lower triangular matrix A, packed columnwise in */
/* > a linear array. The j-th column of A is stored in the array */
/* > AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX*16 array, dimension (N) */
/* > On entry, the right hand side b of the triangular system. */
/* > On exit, X is overwritten by the solution vector x. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* > SCALE is DOUBLE PRECISION */
/* > The scaling factor s for the triangular system */
/* > A * x = s*b, A**T * x = s*b, or A**H * x = s*b. */
/* > If SCALE = 0, the matrix A is singular or badly scaled, and */
/* > the vector x is an exact or approximate solution to A*x = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CNORM */
/* > \verbatim */
/* > CNORM is DOUBLE PRECISION array, dimension (N) */
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
/* > \ingroup latps */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > A rough bound on x is computed;
if that is less than overflow, ZTPSV */
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
/* > Since |x(j)| <= M(j), we use the Level 2 BLAS routine ZTPSV if the */
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
/* > and we can safely call ZTPSV if 1/M(n) and 1/G(n) are both greater */
/* > than fla_max(underflow, 1/overflow). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void zlatps_(char *uplo, char *trans, char *diag, char *normin, integer *n, doublecomplex *ap,
             doublecomplex *x, doublereal *scale, doublereal *cnorm, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlatps inputs: uplo %c, trans %c, diag %c, normin %c, n %" FLA_IS "", *uplo,
                      *trans, *diag, *normin, *n);
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    integer i__, j, ip;
    doublereal xj, rec, tjj;
    integer jinc, jlen;
    doublereal xbnd;
    integer imax;
    doublereal tmax;
    doublecomplex tjjs;
    doublereal xmax, grow;
    extern /* Subroutine */
        void
        dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    doublereal tscal;
    doublecomplex uscal;
    integer jlast;
    doublecomplex csumj;
    extern /* Double Complex */
        void
        zdotc_f2c_(doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                   integer *);
    logical upper;
    extern /* Double Complex */
        void
        zdotu_f2c_(doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                   integer *);
    extern /* Subroutine */
        void
        zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *),
        ztpsv_(char *, char *, char *, integer *, doublecomplex *, doublecomplex *, integer *);
    extern doublereal dlamch_(char *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len),
        zdscal_(integer *, doublereal *, doublecomplex *, integer *);
    doublereal bignum;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    extern /* Double Complex */
        void
        zladiv_f2c_(doublecomplex *, doublecomplex *, doublecomplex *);
    logical notran;
    integer jfirst;
    extern doublereal dzasum_(integer *, doublecomplex *, integer *);
    doublereal smlnum;
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
    --cnorm;
    --x;
    --ap;
    /* Function Body */
    *info = 0;
    // initializing as {1, 0} because it is
    // used as divisor
    tjjs = (doublecomplex){.r = 1., .i = 0.};
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
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZLATPS", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Determine machine dependent parameters to control overflow. */
    smlnum = dlamch_("Safe minimum");
    bignum = 1. / smlnum;
    smlnum /= dlamch_("Precision");
    bignum = 1. / smlnum;
    *scale = 1.;
    if(lsame_(normin, "N", 1, 1))
    {
        /* Compute the 1-norm of each column, not including the diagonal. */
        if(upper)
        {
            /* A is upper triangular. */
            ip = 1;
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = j - 1;
                cnorm[j] = dzasum_(&i__2, &ap[ip], &c__1);
                ip += j;
                /* L10: */
            }
        }
        else
        {
            /* A is lower triangular. */
            ip = 1;
            i__1 = *n - 1;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = *n - j;
                cnorm[j] = dzasum_(&i__2, &ap[ip + 1], &c__1);
                ip = ip + *n - j + 1;
                /* L20: */
            }
            cnorm[*n] = 0.;
        }
    }
    /* Scale the column norms by TSCAL if the maximum element in CNORM is */
    /* greater than BIGNUM/2. */
    imax = idamax_(n, &cnorm[1], &c__1);
    tmax = cnorm[imax];
    if(tmax <= bignum * .5)
    {
        tscal = 1.;
    }
    else
    {
        tscal = .5 / (smlnum * tmax);
        dscal_(n, &tscal, &cnorm[1], &c__1);
    }
    /* Compute a bound on the computed solution vector to see if the */
    /* Level 2 BLAS routine ZTPSV can be used. */
    xmax = 0.;
    i__1 = *n;
    for(j = 1; j <= i__1; ++j)
    {
        /* Computing MAX */
        i__2 = j;
        d__3 = xmax;
        d__4 = (d__1 = x[i__2].r / 2., f2c_abs(d__1))
               + (d__2 = d_imag(&x[j]) / 2., f2c_abs(d__2)); // , expr subst
        xmax = fla_max(d__3, d__4);
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
        if(tscal != 1.)
        {
            grow = 0.;
            goto L60;
        }
        if(nounit)
        {
            /* A is non-unit triangular. */
            /* Compute GROW = 1/G(j) and XBND = 1/M(j). */
            /* Initially, G(0) = max{
           x(i), i=1,...,n}
           . */
            grow = .5 / fla_max(xbnd, smlnum);
            xbnd = grow;
            ip = jfirst * (jfirst + 1) / 2;
            jlen = *n;
            i__1 = jlast;
            i__2 = jinc;
            for(j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
            {
                /* Exit the loop if the growth factor is too small. */
                if(grow <= smlnum)
                {
                    goto L60;
                }
                i__3 = ip;
                tjjs.r = ap[i__3].r;
                tjjs.i = ap[i__3].i; // , expr subst
                tjj = (d__1 = tjjs.r, f2c_abs(d__1)) + (d__2 = d_imag(&tjjs), f2c_abs(d__2));
                if(tjj >= smlnum)
                {
                    /* M(j) = G(j-1) / f2c_abs(A(j,j)) */
                    /* Computing MIN */
                    d__1 = xbnd;
                    d__2 = fla_min(1., tjj) * grow; // , expr subst
                    xbnd = fla_min(d__1, d__2);
                }
                else
                {
                    /* M(j) could overflow, set XBND to 0. */
                    xbnd = 0.;
                }
                if(tjj + cnorm[j] >= smlnum)
                {
                    /* G(j) = G(j-1)*( 1 + CNORM(j) / f2c_abs(A(j,j)) ) */
                    grow *= tjj / (tjj + cnorm[j]);
                }
                else
                {
                    /* G(j) could overflow, set GROW to 0. */
                    grow = 0.;
                }
                ip += jinc * jlen;
                --jlen;
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
            d__1 = 1.;
            d__2 = .5 / fla_max(xbnd, smlnum); // , expr subst
            grow = fla_min(d__1, d__2);
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
                grow *= 1. / (cnorm[j] + 1.);
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
        if(tscal != 1.)
        {
            grow = 0.;
            goto L90;
        }
        if(nounit)
        {
            /* A is non-unit triangular. */
            /* Compute GROW = 1/G(j) and XBND = 1/M(j). */
            /* Initially, M(0) = max{
           x(i), i=1,...,n}
           . */
            grow = .5 / fla_max(xbnd, smlnum);
            xbnd = grow;
            ip = jfirst * (jfirst + 1) / 2;
            jlen = 1;
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
                xj = cnorm[j] + 1.;
                /* Computing MIN */
                d__1 = grow;
                d__2 = xbnd / xj; // , expr subst
                grow = fla_min(d__1, d__2);
                i__3 = ip;
                tjjs.r = ap[i__3].r;
                tjjs.i = ap[i__3].i; // , expr subst
                tjj = (d__1 = tjjs.r, f2c_abs(d__1)) + (d__2 = d_imag(&tjjs), f2c_abs(d__2));
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
                    xbnd = 0.;
                }
                ++jlen;
                ip += jinc * jlen;
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
            d__1 = 1.;
            d__2 = .5 / fla_max(xbnd, smlnum); // , expr subst
            grow = fla_min(d__1, d__2);
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
                xj = cnorm[j] + 1.;
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
        ztpsv_(uplo, trans, diag, n, &ap[1], &x[1], &c__1);
    }
    else
    {
        /* Use a Level 1 BLAS solve, scaling intermediate results. */
        if(xmax > bignum * .5)
        {
            /* Scale X so that its components are less than or equal to */
            /* BIGNUM in absolute value. */
            *scale = bignum * .5 / xmax;
            zdscal_(n, scale, &x[1], &c__1);
            xmax = bignum;
        }
        else
        {
            xmax *= 2.;
        }
        if(notran)
        {
            /* Solve A * x = b */
            ip = jfirst * (jfirst + 1) / 2;
            i__1 = jlast;
            i__2 = jinc;
            for(j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
            {
                /* Compute x(j) = b(j) / A(j,j), scaling x if necessary. */
                i__3 = j;
                xj = (d__1 = x[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&x[j]), f2c_abs(d__2));
                if(nounit)
                {
                    i__3 = ip;
                    z__1.r = tscal * ap[i__3].r;
                    z__1.i = tscal * ap[i__3].i; // , expr subst
                    tjjs.r = z__1.r;
                    tjjs.i = z__1.i; // , expr subst
                }
                else
                {
                    tjjs.r = tscal;
                    tjjs.i = 0.; // , expr subst
                    if(tscal == 1.)
                    {
                        goto L110;
                    }
                }
                tjj = (d__1 = tjjs.r, f2c_abs(d__1)) + (d__2 = d_imag(&tjjs), f2c_abs(d__2));
                if(tjj > smlnum)
                {
                    /* f2c_abs(A(j,j)) > SMLNUM: */
                    if(tjj < 1.)
                    {
                        if(xj > tjj * bignum)
                        {
                            /* Scale x by 1/b(j). */
                            rec = 1. / xj;
                            zdscal_(n, &rec, &x[1], &c__1);
                            *scale *= rec;
                            xmax *= rec;
                        }
                    }
                    i__3 = j;
                    zladiv_f2c_(&z__1, &x[j], &tjjs);
                    x[i__3].r = z__1.r;
                    x[i__3].i = z__1.i; // , expr subst
                    i__3 = j;
                    xj = (d__1 = x[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&x[j]), f2c_abs(d__2));
                }
                else if(tjj > 0.)
                {
                    /* 0 < f2c_abs(A(j,j)) <= SMLNUM: */
                    if(xj > tjj * bignum)
                    {
                        /* Scale x by (1/f2c_abs(x(j)))*f2c_abs(A(j,j))*BIGNUM */
                        /* to avoid overflow when dividing by A(j,j). */
                        rec = tjj * bignum / xj;
                        if(cnorm[j] > 1.)
                        {
                            /* Scale by 1/CNORM(j) to avoid overflow when */
                            /* multiplying x(j) times column j. */
                            rec /= cnorm[j];
                        }
                        zdscal_(n, &rec, &x[1], &c__1);
                        *scale *= rec;
                        xmax *= rec;
                    }
                    i__3 = j;
                    zladiv_f2c_(&z__1, &x[j], &tjjs);
                    x[i__3].r = z__1.r;
                    x[i__3].i = z__1.i; // , expr subst
                    i__3 = j;
                    xj = (d__1 = x[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&x[j]), f2c_abs(d__2));
                }
                else
                {
                    /* A(j,j) = 0: Set x(1:n) = 0, x(j) = 1, and */
                    /* scale = 0, and compute a solution to A*x = 0. */
                    i__3 = *n;
                    for(i__ = 1; i__ <= i__3; ++i__)
                    {
                        i__4 = i__;
                        x[i__4].r = 0.;
                        x[i__4].i = 0.; // , expr subst
                        /* L100: */
                    }
                    i__3 = j;
                    x[i__3].r = 1.;
                    x[i__3].i = 0.; // , expr subst
                    xj = 1.;
                    *scale = 0.;
                    xmax = 0.;
                }
            L110: /* Scale x if necessary to avoid overflow when adding a */
                /* multiple of column j of A. */
                if(xj > 1.)
                {
                    rec = 1. / xj;
                    if(cnorm[j] > (bignum - xmax) * rec)
                    {
                        /* Scale x by 1/(2*f2c_abs(x(j))). */
                        rec *= .5;
                        zdscal_(n, &rec, &x[1], &c__1);
                        *scale *= rec;
                    }
                }
                else if(xj * cnorm[j] > bignum - xmax)
                {
                    /* Scale x by 1/2. */
                    zdscal_(n, &c_b36, &x[1], &c__1);
                    *scale *= .5;
                }
                if(upper)
                {
                    if(j > 1)
                    {
                        /* Compute the update */
                        /* x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j) */
                        i__3 = j - 1;
                        i__4 = j;
                        z__2.r = -x[i__4].r;
                        z__2.i = -x[i__4].i; // , expr subst
                        z__1.r = tscal * z__2.r;
                        z__1.i = tscal * z__2.i; // , expr subst
                        zaxpy_(&i__3, &z__1, &ap[ip - j + 1], &c__1, &x[1], &c__1);
                        i__3 = j - 1;
                        i__ = izamax_(&i__3, &x[1], &c__1);
                        i__3 = i__;
                        xmax = (d__1 = x[i__3].r, f2c_abs(d__1))
                               + (d__2 = d_imag(&x[i__]), f2c_abs(d__2));
                    }
                    ip -= j;
                }
                else
                {
                    if(j < *n)
                    {
                        /* Compute the update */
                        /* x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j) */
                        i__3 = *n - j;
                        i__4 = j;
                        z__2.r = -x[i__4].r;
                        z__2.i = -x[i__4].i; // , expr subst
                        z__1.r = tscal * z__2.r;
                        z__1.i = tscal * z__2.i; // , expr subst
                        zaxpy_(&i__3, &z__1, &ap[ip + 1], &c__1, &x[j + 1], &c__1);
                        i__3 = *n - j;
                        i__ = j + izamax_(&i__3, &x[j + 1], &c__1);
                        i__3 = i__;
                        xmax = (d__1 = x[i__3].r, f2c_abs(d__1))
                               + (d__2 = d_imag(&x[i__]), f2c_abs(d__2));
                    }
                    ip = ip + *n - j + 1;
                }
                /* L120: */
            }
        }
        else if(lsame_(trans, "T", 1, 1))
        {
            /* Solve A**T * x = b */
            ip = jfirst * (jfirst + 1) / 2;
            jlen = 1;
            i__2 = jlast;
            i__1 = jinc;
            for(j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1)
            {
                /* Compute x(j) = b(j) - sum A(k,j)*x(k). */
                /* k<>j */
                i__3 = j;
                xj = (d__1 = x[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&x[j]), f2c_abs(d__2));
                uscal.r = tscal;
                uscal.i = 0.; // , expr subst
                rec = 1. / fla_max(xmax, 1.);
                if(cnorm[j] > (bignum - xj) * rec)
                {
                    /* If x(j) could overflow, scale x by 1/(2*XMAX). */
                    rec *= .5;
                    if(nounit)
                    {
                        i__3 = ip;
                        z__1.r = tscal * ap[i__3].r;
                        z__1.i = tscal * ap[i__3].i; // , expr subst
                        tjjs.r = z__1.r;
                        tjjs.i = z__1.i; // , expr subst
                    }
                    else
                    {
                        tjjs.r = tscal;
                        tjjs.i = 0.; // , expr subst
                    }
                    tjj = (d__1 = tjjs.r, f2c_abs(d__1)) + (d__2 = d_imag(&tjjs), f2c_abs(d__2));
                    if(tjj > 1.)
                    {
                        /* Divide by A(j,j) when scaling x if A(j,j) > 1. */
                        /* Computing MIN */
                        d__1 = 1.;
                        d__2 = rec * tjj; // , expr subst
                        rec = fla_min(d__1, d__2);
                        zladiv_f2c_(&z__1, &uscal, &tjjs);
                        uscal.r = z__1.r;
                        uscal.i = z__1.i; // , expr subst
                    }
                    if(rec < 1.)
                    {
                        zdscal_(n, &rec, &x[1], &c__1);
                        *scale *= rec;
                        xmax *= rec;
                    }
                }
                csumj.r = 0.;
                csumj.i = 0.; // , expr subst
                if(uscal.r == 1. && uscal.i == 0.)
                {
                    /* If the scaling needed for A in the dot product is 1, */
                    /* call ZDOTU to perform the dot product. */
                    if(upper)
                    {
                        i__3 = j - 1;
                        zdotu_f2c_(&z__1, &i__3, &ap[ip - j + 1], &c__1, &x[1], &c__1);
                        csumj.r = z__1.r;
                        csumj.i = z__1.i; // , expr subst
                    }
                    else if(j < *n)
                    {
                        i__3 = *n - j;
                        zdotu_f2c_(&z__1, &i__3, &ap[ip + 1], &c__1, &x[j + 1], &c__1);
                        csumj.r = z__1.r;
                        csumj.i = z__1.i; // , expr subst
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
                            i__4 = ip - j + i__;
                            z__3.r = ap[i__4].r * uscal.r - ap[i__4].i * uscal.i;
                            z__3.i = ap[i__4].r * uscal.i + ap[i__4].i * uscal.r; // , expr subst
                            i__5 = i__;
                            z__2.r = z__3.r * x[i__5].r - z__3.i * x[i__5].i;
                            z__2.i = z__3.r * x[i__5].i + z__3.i * x[i__5].r; // , expr subst
                            z__1.r = csumj.r + z__2.r;
                            z__1.i = csumj.i + z__2.i; // , expr subst
                            csumj.r = z__1.r;
                            csumj.i = z__1.i; // , expr subst
                            /* L130: */
                        }
                    }
                    else if(j < *n)
                    {
                        i__3 = *n - j;
                        for(i__ = 1; i__ <= i__3; ++i__)
                        {
                            i__4 = ip + i__;
                            z__3.r = ap[i__4].r * uscal.r - ap[i__4].i * uscal.i;
                            z__3.i = ap[i__4].r * uscal.i + ap[i__4].i * uscal.r; // , expr subst
                            i__5 = j + i__;
                            z__2.r = z__3.r * x[i__5].r - z__3.i * x[i__5].i;
                            z__2.i = z__3.r * x[i__5].i + z__3.i * x[i__5].r; // , expr subst
                            z__1.r = csumj.r + z__2.r;
                            z__1.i = csumj.i + z__2.i; // , expr subst
                            csumj.r = z__1.r;
                            csumj.i = z__1.i; // , expr subst
                            /* L140: */
                        }
                    }
                }
                z__1.r = tscal;
                z__1.i = 0.; // , expr subst
                if(uscal.r == z__1.r && uscal.i == z__1.i)
                {
                    /* Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j) */
                    /* was not used to scale the dotproduct. */
                    i__3 = j;
                    i__4 = j;
                    z__1.r = x[i__4].r - csumj.r;
                    z__1.i = x[i__4].i - csumj.i; // , expr subst
                    x[i__3].r = z__1.r;
                    x[i__3].i = z__1.i; // , expr subst
                    i__3 = j;
                    xj = (d__1 = x[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&x[j]), f2c_abs(d__2));
                    if(nounit)
                    {
                        /* Compute x(j) = x(j) / A(j,j), scaling if necessary. */
                        i__3 = ip;
                        z__1.r = tscal * ap[i__3].r;
                        z__1.i = tscal * ap[i__3].i; // , expr subst
                        tjjs.r = z__1.r;
                        tjjs.i = z__1.i; // , expr subst
                    }
                    else
                    {
                        tjjs.r = tscal;
                        tjjs.i = 0.; // , expr subst
                        if(tscal == 1.)
                        {
                            goto L160;
                        }
                    }
                    tjj = (d__1 = tjjs.r, f2c_abs(d__1)) + (d__2 = d_imag(&tjjs), f2c_abs(d__2));
                    if(tjj > smlnum)
                    {
                        /* f2c_abs(A(j,j)) > SMLNUM: */
                        if(tjj < 1.)
                        {
                            if(xj > tjj * bignum)
                            {
                                /* Scale X by 1/f2c_abs(x(j)). */
                                rec = 1. / xj;
                                zdscal_(n, &rec, &x[1], &c__1);
                                *scale *= rec;
                                xmax *= rec;
                            }
                        }
                        i__3 = j;
                        zladiv_f2c_(&z__1, &x[j], &tjjs);
                        x[i__3].r = z__1.r;
                        x[i__3].i = z__1.i; // , expr subst
                    }
                    else if(tjj > 0.)
                    {
                        /* 0 < f2c_abs(A(j,j)) <= SMLNUM: */
                        if(xj > tjj * bignum)
                        {
                            /* Scale x by (1/f2c_abs(x(j)))*f2c_abs(A(j,j))*BIGNUM. */
                            rec = tjj * bignum / xj;
                            zdscal_(n, &rec, &x[1], &c__1);
                            *scale *= rec;
                            xmax *= rec;
                        }
                        i__3 = j;
                        zladiv_f2c_(&z__1, &x[j], &tjjs);
                        x[i__3].r = z__1.r;
                        x[i__3].i = z__1.i; // , expr subst
                    }
                    else
                    {
                        /* A(j,j) = 0: Set x(1:n) = 0, x(j) = 1, and */
                        /* scale = 0 and compute a solution to A**T *x = 0. */
                        i__3 = *n;
                        for(i__ = 1; i__ <= i__3; ++i__)
                        {
                            i__4 = i__;
                            x[i__4].r = 0.;
                            x[i__4].i = 0.; // , expr subst
                            /* L150: */
                        }
                        i__3 = j;
                        x[i__3].r = 1.;
                        x[i__3].i = 0.; // , expr subst
                        *scale = 0.;
                        xmax = 0.;
                    }
                L160:;
                }
                else
                {
                    /* Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot */
                    /* product has already been divided by 1/A(j,j). */
                    i__3 = j;
                    zladiv_f2c_(&z__2, &x[j], &tjjs);
                    z__1.r = z__2.r - csumj.r;
                    z__1.i = z__2.i - csumj.i; // , expr subst
                    x[i__3].r = z__1.r;
                    x[i__3].i = z__1.i; // , expr subst
                }
                /* Computing MAX */
                i__3 = j;
                d__3 = xmax;
                d__4 = (d__1 = x[i__3].r, f2c_abs(d__1))
                       + (d__2 = d_imag(&x[j]), f2c_abs(d__2)); // , expr subst
                xmax = fla_max(d__3, d__4);
                ++jlen;
                ip += jinc * jlen;
                /* L170: */
            }
        }
        else
        {
            /* Solve A**H * x = b */
            ip = jfirst * (jfirst + 1) / 2;
            jlen = 1;
            i__1 = jlast;
            i__2 = jinc;
            for(j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
            {
                /* Compute x(j) = b(j) - sum A(k,j)*x(k). */
                /* k<>j */
                i__3 = j;
                xj = (d__1 = x[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&x[j]), f2c_abs(d__2));
                uscal.r = tscal;
                uscal.i = 0.; // , expr subst
                rec = 1. / fla_max(xmax, 1.);
                if(cnorm[j] > (bignum - xj) * rec)
                {
                    /* If x(j) could overflow, scale x by 1/(2*XMAX). */
                    rec *= .5;
                    if(nounit)
                    {
                        d_cnjg(&z__2, &ap[ip]);
                        z__1.r = tscal * z__2.r;
                        z__1.i = tscal * z__2.i; // , expr subst
                        tjjs.r = z__1.r;
                        tjjs.i = z__1.i; // , expr subst
                    }
                    else
                    {
                        tjjs.r = tscal;
                        tjjs.i = 0.; // , expr subst
                    }
                    tjj = (d__1 = tjjs.r, f2c_abs(d__1)) + (d__2 = d_imag(&tjjs), f2c_abs(d__2));
                    if(tjj > 1.)
                    {
                        /* Divide by A(j,j) when scaling x if A(j,j) > 1. */
                        /* Computing MIN */
                        d__1 = 1.;
                        d__2 = rec * tjj; // , expr subst
                        rec = fla_min(d__1, d__2);
                        zladiv_f2c_(&z__1, &uscal, &tjjs);
                        uscal.r = z__1.r;
                        uscal.i = z__1.i; // , expr subst
                    }
                    if(rec < 1.)
                    {
                        zdscal_(n, &rec, &x[1], &c__1);
                        *scale *= rec;
                        xmax *= rec;
                    }
                }
                csumj.r = 0.;
                csumj.i = 0.; // , expr subst
                if(uscal.r == 1. && uscal.i == 0.)
                {
                    /* If the scaling needed for A in the dot product is 1, */
                    /* call ZDOTC to perform the dot product. */
                    if(upper)
                    {
                        i__3 = j - 1;
                        zdotc_f2c_(&z__1, &i__3, &ap[ip - j + 1], &c__1, &x[1], &c__1);
                        csumj.r = z__1.r;
                        csumj.i = z__1.i; // , expr subst
                    }
                    else if(j < *n)
                    {
                        i__3 = *n - j;
                        zdotc_f2c_(&z__1, &i__3, &ap[ip + 1], &c__1, &x[j + 1], &c__1);
                        csumj.r = z__1.r;
                        csumj.i = z__1.i; // , expr subst
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
                            d_cnjg(&z__4, &ap[ip - j + i__]);
                            z__3.r = z__4.r * uscal.r - z__4.i * uscal.i;
                            z__3.i = z__4.r * uscal.i + z__4.i * uscal.r; // , expr subst
                            i__4 = i__;
                            z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i;
                            z__2.i = z__3.r * x[i__4].i + z__3.i * x[i__4].r; // , expr subst
                            z__1.r = csumj.r + z__2.r;
                            z__1.i = csumj.i + z__2.i; // , expr subst
                            csumj.r = z__1.r;
                            csumj.i = z__1.i; // , expr subst
                            /* L180: */
                        }
                    }
                    else if(j < *n)
                    {
                        i__3 = *n - j;
                        for(i__ = 1; i__ <= i__3; ++i__)
                        {
                            d_cnjg(&z__4, &ap[ip + i__]);
                            z__3.r = z__4.r * uscal.r - z__4.i * uscal.i;
                            z__3.i = z__4.r * uscal.i + z__4.i * uscal.r; // , expr subst
                            i__4 = j + i__;
                            z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i;
                            z__2.i = z__3.r * x[i__4].i + z__3.i * x[i__4].r; // , expr subst
                            z__1.r = csumj.r + z__2.r;
                            z__1.i = csumj.i + z__2.i; // , expr subst
                            csumj.r = z__1.r;
                            csumj.i = z__1.i; // , expr subst
                            /* L190: */
                        }
                    }
                }
                z__1.r = tscal;
                z__1.i = 0.; // , expr subst
                if(uscal.r == z__1.r && uscal.i == z__1.i)
                {
                    /* Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j) */
                    /* was not used to scale the dotproduct. */
                    i__3 = j;
                    i__4 = j;
                    z__1.r = x[i__4].r - csumj.r;
                    z__1.i = x[i__4].i - csumj.i; // , expr subst
                    x[i__3].r = z__1.r;
                    x[i__3].i = z__1.i; // , expr subst
                    i__3 = j;
                    xj = (d__1 = x[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&x[j]), f2c_abs(d__2));
                    if(nounit)
                    {
                        /* Compute x(j) = x(j) / A(j,j), scaling if necessary. */
                        d_cnjg(&z__2, &ap[ip]);
                        z__1.r = tscal * z__2.r;
                        z__1.i = tscal * z__2.i; // , expr subst
                        tjjs.r = z__1.r;
                        tjjs.i = z__1.i; // , expr subst
                    }
                    else
                    {
                        tjjs.r = tscal;
                        tjjs.i = 0.; // , expr subst
                        if(tscal == 1.)
                        {
                            goto L210;
                        }
                    }
                    tjj = (d__1 = tjjs.r, f2c_abs(d__1)) + (d__2 = d_imag(&tjjs), f2c_abs(d__2));
                    if(tjj > smlnum)
                    {
                        /* f2c_abs(A(j,j)) > SMLNUM: */
                        if(tjj < 1.)
                        {
                            if(xj > tjj * bignum)
                            {
                                /* Scale X by 1/f2c_abs(x(j)). */
                                rec = 1. / xj;
                                zdscal_(n, &rec, &x[1], &c__1);
                                *scale *= rec;
                                xmax *= rec;
                            }
                        }
                        i__3 = j;
                        zladiv_f2c_(&z__1, &x[j], &tjjs);
                        x[i__3].r = z__1.r;
                        x[i__3].i = z__1.i; // , expr subst
                    }
                    else if(tjj > 0.)
                    {
                        /* 0 < f2c_abs(A(j,j)) <= SMLNUM: */
                        if(xj > tjj * bignum)
                        {
                            /* Scale x by (1/f2c_abs(x(j)))*f2c_abs(A(j,j))*BIGNUM. */
                            rec = tjj * bignum / xj;
                            zdscal_(n, &rec, &x[1], &c__1);
                            *scale *= rec;
                            xmax *= rec;
                        }
                        i__3 = j;
                        zladiv_f2c_(&z__1, &x[j], &tjjs);
                        x[i__3].r = z__1.r;
                        x[i__3].i = z__1.i; // , expr subst
                    }
                    else
                    {
                        /* A(j,j) = 0: Set x(1:n) = 0, x(j) = 1, and */
                        /* scale = 0 and compute a solution to A**H *x = 0. */
                        i__3 = *n;
                        for(i__ = 1; i__ <= i__3; ++i__)
                        {
                            i__4 = i__;
                            x[i__4].r = 0.;
                            x[i__4].i = 0.; // , expr subst
                            /* L200: */
                        }
                        i__3 = j;
                        x[i__3].r = 1.;
                        x[i__3].i = 0.; // , expr subst
                        *scale = 0.;
                        xmax = 0.;
                    }
                L210:;
                }
                else
                {
                    /* Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot */
                    /* product has already been divided by 1/A(j,j). */
                    i__3 = j;
                    zladiv_f2c_(&z__2, &x[j], &tjjs);
                    z__1.r = z__2.r - csumj.r;
                    z__1.i = z__2.i - csumj.i; // , expr subst
                    x[i__3].r = z__1.r;
                    x[i__3].i = z__1.i; // , expr subst
                }
                /* Computing MAX */
                i__3 = j;
                d__3 = xmax;
                d__4 = (d__1 = x[i__3].r, f2c_abs(d__1))
                       + (d__2 = d_imag(&x[j]), f2c_abs(d__2)); // , expr subst
                xmax = fla_max(d__3, d__4);
                ++jlen;
                ip += jinc * jlen;
                /* L220: */
            }
        }
        *scale /= tscal;
    }
    /* Scale the column norms by 1/TSCAL for return. */
    if(tscal != 1.)
    {
        d__1 = 1. / tscal;
        dscal_(n, &d__1, &cnorm[1], &c__1);
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZLATPS */
}
/* zlatps_ */
