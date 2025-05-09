/* ./zbdsqr.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b15 = -.125;
static integer c__1 = 1;
static doublereal c_b49 = 1.;
static doublereal c_b72 = -1.;
/* > \brief \b ZBDSQR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZBDSQR + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zbdsqr.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zbdsqr.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zbdsqr.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, */
/* LDU, C, LDC, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION D( * ), E( * ), RWORK( * ) */
/* COMPLEX*16 C( LDC, * ), U( LDU, * ), VT( LDVT, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZBDSQR computes the singular values and, optionally, the right and/or */
/* > left singular vectors from the singular value decomposition (SVD) of */
/* > a real N-by-N (upper or lower) bidiagonal matrix B using the implicit */
/* > zero-shift QR algorithm. The SVD of B has the form */
/* > */
/* > B = Q * S * P**H */
/* > */
/* > where S is the diagonal matrix of singular values, Q is an orthogonal */
/* > matrix of left singular vectors, and P is an orthogonal matrix of */
/* > right singular vectors. If left singular vectors are requested, this */
/* > subroutine actually returns U*Q instead of Q, and, if right singular */
/* > vectors are requested, this subroutine returns P**H*VT instead of */
/* > P**H, for given complex input matrices U and VT. When U and VT are */
/* > the unitary matrices that reduce a general matrix A to bidiagonal */
/* > form: A = U*B*VT, as computed by ZGEBRD, then */
/* > */
/* > A = (U*Q) * S * (P**H*VT) */
/* > */
/* > is the SVD of A. Optionally, the subroutine may also compute Q**H*C */
/* > for a given complex input matrix C. */
/* > */
/* > See "Computing Small Singular Values of Bidiagonal Matrices With */
/* > Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan, */
/* > LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11, */
/* > no. 5, pp. 873-912, Sept 1990) and */
/* > "Accurate singular values and differential qd algorithms," by */
/* > B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics */
/* > Department, University of California at Berkeley, July 1992 */
/* > for a detailed description of the algorithm. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': B is upper bidiagonal;
 */
/* > = 'L': B is lower bidiagonal. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NCVT */
/* > \verbatim */
/* > NCVT is INTEGER */
/* > The number of columns of the matrix VT. NCVT >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRU */
/* > \verbatim */
/* > NRU is INTEGER */
/* > The number of rows of the matrix U. NRU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NCC */
/* > \verbatim */
/* > NCC is INTEGER */
/* > The number of columns of the matrix C. NCC >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (N) */
/* > On entry, the n diagonal elements of the bidiagonal matrix B. */
/* > On exit, if INFO=0, the singular values of B in decreasing */
/* > order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is DOUBLE PRECISION array, dimension (N-1) */
/* > On entry, the N-1 offdiagonal elements of the bidiagonal */
/* > matrix B. */
/* > On exit, if INFO = 0, E is destroyed;
if INFO > 0, D and E */
/* > will contain the diagonal and superdiagonal elements of a */
/* > bidiagonal matrix orthogonally equivalent to the one given */
/* > as input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VT */
/* > \verbatim */
/* > VT is COMPLEX*16 array, dimension (LDVT, NCVT) */
/* > On entry, an N-by-NCVT matrix VT. */
/* > On exit, VT is overwritten by P**H * VT. */
/* > Not referenced if NCVT = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* > LDVT is INTEGER */
/* > The leading dimension of the array VT. */
/* > LDVT >= fla_max(1,N) if NCVT > 0;
LDVT >= 1 if NCVT = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] U */
/* > \verbatim */
/* > U is COMPLEX*16 array, dimension (LDU, N) */
/* > On entry, an NRU-by-N matrix U. */
/* > On exit, U is overwritten by U * Q. */
/* > Not referenced if NRU = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* > LDU is INTEGER */
/* > The leading dimension of the array U. LDU >= fla_max(1,NRU). */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX*16 array, dimension (LDC, NCC) */
/* > On entry, an N-by-NCC matrix C. */
/* > On exit, C is overwritten by Q**H * C. */
/* > Not referenced if NCC = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. */
/* > LDC >= fla_max(1,N) if NCC > 0;
LDC >=1 if NCC = 0. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is DOUBLE PRECISION array, dimension (4*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: If INFO = -i, the i-th argument had an illegal value */
/* > > 0: the algorithm did not converge;
D and E contain the */
/* > elements of a bidiagonal matrix which is orthogonally */
/* > similar to the input matrix B;
if INFO = i, i */
/* > elements of E have not converged to zero. */
/* > \endverbatim */
/* > \par Internal Parameters: */
/* ========================= */
/* > */
/* > \verbatim */
/* > TOLMUL DOUBLE PRECISION, default = fla_max(10,min(100,EPS**(-1/8))) */
/* > TOLMUL controls the convergence criterion of the QR loop. */
/* > If it is positive, TOLMUL*EPS is the desired relative */
/* > precision in the computed singular values. */
/* > If it is negative, f2c_dabs(TOLMUL*EPS*sigma_max) is the */
/* > desired absolute accuracy in the computed singular */
/* > values (corresponds to relative accuracy */
/* > f2c_dabs(TOLMUL*EPS) in the largest singular value. */
/* > f2c_dabs(TOLMUL) should be between 1 and 1/EPS, and preferably */
/* > between 10 (for fast convergence) and .1/EPS */
/* > (for there to be some accuracy in the results). */
/* > Default is to lose at either one eighth or 2 of the */
/* > available decimal digits in each computed singular value */
/* > (whichever is smaller). */
/* > */
/* > MAXITR INTEGER, default = 6 */
/* > MAXITR controls the maximum number of passes of the */
/* > algorithm through its inner loop. The algorithms stops */
/* > (and so fails to converge) if the number of passes */
/* > through the inner loop exceeds MAXITR*N**2. */
/* > */
/* > \endverbatim */
/* > \par Note: */
/* =========== */
/* > */
/* > \verbatim */
/* > Bug report from Cezary Dendek. */
/* > On November 3rd 2023, the INTEGER variable MAXIT = MAXITR*N**2 is */
/* > removed since it can overflow pretty easily (for N larger or equal */
/* > than 18,919). We instead use MAXITDIVN = MAXITR*N. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup bdsqr */
/* ===================================================================== */
/* Subroutine */
void zbdsqr_(char *uplo, integer *n, integer *ncvt, integer *nru, integer *ncc, doublereal *d__,
             doublereal *e, doublecomplex *vt, integer *ldvt, doublecomplex *u, integer *ldu,
             doublecomplex *c__, integer *ldc, doublereal *rwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zbdsqr inputs: uplo %c, n %" FLA_IS ", ncvt %" FLA_IS ", nru %" FLA_IS
                      ", ncc %" FLA_IS ", ldvt %" FLA_IS ", ldu %" FLA_IS ", ldc %" FLA_IS "",
                      *uplo, *n, *ncvt, *nru, *ncc, *ldvt, *ldu, *ldc);
    /* System generated locals */
    integer c_dim1, c_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;
    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), d_sign(doublereal *, doublereal *);
    /* Local variables */
    integer iterdivn;
    doublereal f, g, h__;
    integer i__, j, m;
    doublereal r__;
    integer maxitdivn;
    doublereal cs;
    integer ll;
    doublereal sn, mu;
    integer nm1, nm12, nm13, lll;
    doublereal eps, sll, tol, abse;
    integer idir;
    doublereal abss;
    integer oldm;
    doublereal cosl;
    integer isub, iter;
    doublereal unfl, sinl, cosr, smin, smax, sinr;
    extern /* Subroutine */
        void
        dlas2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    extern logical lsame_(char *, char *, integer, integer);
    doublereal oldcs;
    integer oldll;
    doublereal shift, sigmn, oldsn, sigmx;
    logical lower;
    extern /* Subroutine */
        void
        zlasr_(char *, char *, char *, integer *, integer *, doublereal *, doublereal *,
               doublecomplex *, integer *),
        zdrot_(integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublereal *,
               doublereal *),
        zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *),
        dlasq1_(integer *, doublereal *, doublereal *, doublereal *, integer *),
        dlasv2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, doublereal *, doublereal *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
        void
        dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len),
        zdscal_(integer *, doublereal *, doublecomplex *, integer *);
    doublereal sminoa, thresh;
    logical rotate;
    doublereal tolmul;
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --d__;
    --e;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --rwork;
    /* Function Body */
    *info = 0;
    lower = lsame_(uplo, "L", 1, 1);
    oldsn = 0.;
    if(!lsame_(uplo, "U", 1, 1) && !lower)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*ncvt < 0)
    {
        *info = -3;
    }
    else if(*nru < 0)
    {
        *info = -4;
    }
    else if(*ncc < 0)
    {
        *info = -5;
    }
    else if(*ncvt == 0 && *ldvt < 1 || *ncvt > 0 && *ldvt < fla_max(1, *n))
    {
        *info = -9;
    }
    else if(*ldu < fla_max(1, *nru))
    {
        *info = -11;
    }
    else if(*ncc == 0 && *ldc < 1 || *ncc > 0 && *ldc < fla_max(1, *n))
    {
        *info = -13;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZBDSQR", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*n == 1)
    {
        goto L160;
    }
    /* ROTATE is true if any singular vectors desired, false otherwise */
    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;
    /* If no singular vectors desired, use qd algorithm */
    if(!rotate)
    {
        dlasq1_(n, &d__[1], &e[1], &rwork[1], info);
        /* If INFO equals 2, dqds didn't finish, try to finish */
        if(*info != 2)
        {
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        *info = 0;
    }
    nm1 = *n - 1;
    nm12 = nm1 + nm1;
    nm13 = nm12 + nm1;
    idir = 0;
    /* Get machine constants */
    eps = dlamch_("Epsilon");
    unfl = dlamch_("Safe minimum");
    /* If matrix lower bidiagonal, rotate to be upper bidiagonal */
    /* by applying Givens rotations on the left */
    if(lower)
    {
        i__1 = *n - 1;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
            d__[i__] = r__;
            e[i__] = sn * d__[i__ + 1];
            d__[i__ + 1] = cs * d__[i__ + 1];
            rwork[i__] = cs;
            rwork[nm1 + i__] = sn;
            /* L10: */
        }
        /* Update singular vectors if desired */
        if(*nru > 0)
        {
            zlasr_("R", "V", "F", nru, n, &rwork[1], &rwork[*n], &u[u_offset], ldu);
        }
        if(*ncc > 0)
        {
            zlasr_("L", "V", "F", n, ncc, &rwork[1], &rwork[*n], &c__[c_offset], ldc);
        }
    }
    /* Compute singular values to relative accuracy TOL */
    /* (By setting TOL to be negative, algorithm will compute */
    /* singular values to absolute accuracy ABS(TOL)*norm(input matrix)) */
    /* Computing MAX */
    /* Computing MIN */
    d__3 = 100.;
    d__4 = pow_dd(&eps, &c_b15); // , expr subst
    d__1 = 10.;
    d__2 = fla_min(d__3, d__4); // , expr subst
    tolmul = fla_max(d__1, d__2);
    tol = tolmul * eps;
    /* Compute approximate maximum, minimum singular values */
    smax = 0.;
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        /* Computing MAX */
        d__2 = smax;
        d__3 = (d__1 = d__[i__], f2c_dabs(d__1)); // , expr subst
        smax = fla_max(d__2, d__3);
        /* L20: */
    }
    i__1 = *n - 1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        /* Computing MAX */
        d__2 = smax;
        d__3 = (d__1 = e[i__], f2c_dabs(d__1)); // , expr subst
        smax = fla_max(d__2, d__3);
        /* L30: */
    }
    smin = 0.;
    if(tol >= 0.)
    {
        /* Relative accuracy desired */
        sminoa = f2c_dabs(d__[1]);
        if(sminoa == 0.)
        {
            goto L50;
        }
        mu = sminoa;
        i__1 = *n;
        for(i__ = 2; i__ <= i__1; ++i__)
        {
            mu = (d__2 = d__[i__], f2c_dabs(d__2))
                 * (mu / (mu + (d__1 = e[i__ - 1], f2c_dabs(d__1))));
            sminoa = fla_min(sminoa, mu);
            if(sminoa == 0.)
            {
                goto L50;
            }
            /* L40: */
        }
    L50:
        sminoa /= sqrt((doublereal)(*n));
        /* Computing MAX */
        d__1 = tol * sminoa;
        d__2 = *n * (*n * unfl) * 6; // , expr subst
        thresh = fla_max(d__1, d__2);
    }
    else
    {
        /* Absolute accuracy desired */
        /* Computing MAX */
        d__1 = f2c_dabs(tol) * smax;
        d__2 = *n * (*n * unfl) * 6; // , expr subst
        thresh = fla_max(d__1, d__2);
    }
    /* Prepare for main iteration loop for the singular values */
    /* (MAXIT is the maximum number of passes through the inner */
    /* loop permitted before nonconvergence signalled.) */
    maxitdivn = *n * 6;
    iterdivn = 0;
    iter = -1;
    oldll = -1;
    oldm = -1;
    /* M points to last element of unconverged part of matrix */
    m = *n;
/* Begin main iteration loop */
L60: /* Check for convergence or exceeding iteration count */
    if(m <= 1)
    {
        goto L160;
    }
    if(iter >= *n)
    {
        iter -= *n;
        ++iterdivn;
        if(iterdivn >= maxitdivn)
        {
            goto L200;
        }
    }
    /* Find diagonal block of matrix to work on */
    if(tol < 0. && (d__1 = d__[m], f2c_dabs(d__1)) <= thresh)
    {
        d__[m] = 0.;
    }
    smax = (d__1 = d__[m], f2c_dabs(d__1));
    i__1 = m - 1;
    for(lll = 1; lll <= i__1; ++lll)
    {
        ll = m - lll;
        abss = (d__1 = d__[ll], f2c_dabs(d__1));
        abse = (d__1 = e[ll], f2c_dabs(d__1));
        if(tol < 0. && abss <= thresh)
        {
            d__[ll] = 0.;
        }
        if(abse <= thresh)
        {
            goto L80;
        }
        /* Computing MAX */
        d__1 = fla_max(smax, abss);
        smax = fla_max(d__1, abse);
        /* L70: */
    }
    ll = 0;
    goto L90;
L80:
    e[ll] = 0.;
    /* Matrix splits since E(LL) = 0 */
    if(ll == m - 1)
    {
        /* Convergence of bottom singular value, return to top of loop */
        --m;
        goto L60;
    }
L90:
    ++ll;
    /* E(LL) through E(M-1) are nonzero, E(LL-1) is zero */
    if(ll == m - 1)
    {
        /* 2 by 2 block, handle separately */
        dlasv2_(&d__[m - 1], &e[m - 1], &d__[m], &sigmn, &sigmx, &sinr, &cosr, &sinl, &cosl);
        d__[m - 1] = sigmx;
        e[m - 1] = 0.;
        d__[m] = sigmn;
        /* Compute singular vectors, if desired */
        if(*ncvt > 0)
        {
            zdrot_(ncvt, &vt[m - 1 + vt_dim1], ldvt, &vt[m + vt_dim1], ldvt, &cosr, &sinr);
        }
        if(*nru > 0)
        {
            zdrot_(nru, &u[(m - 1) * u_dim1 + 1], &c__1, &u[m * u_dim1 + 1], &c__1, &cosl, &sinl);
        }
        if(*ncc > 0)
        {
            zdrot_(ncc, &c__[m - 1 + c_dim1], ldc, &c__[m + c_dim1], ldc, &cosl, &sinl);
        }
        m += -2;
        goto L60;
    }
    /* If working on new submatrix, choose shift direction */
    /* (from larger end diagonal element towards smaller) */
    if(ll > oldm || m < oldll)
    {
        if((d__1 = d__[ll], f2c_dabs(d__1)) >= (d__2 = d__[m], f2c_dabs(d__2)))
        {
            /* Chase bulge from top (big end) to bottom (small end) */
            idir = 1;
        }
        else
        {
            /* Chase bulge from bottom (big end) to top (small end) */
            idir = 2;
        }
    }
    /* Apply convergence tests */
    if(idir == 1)
    {
        /* Run convergence test in forward direction */
        /* First apply standard test to bottom of matrix */
        if((d__2 = e[m - 1], f2c_dabs(d__2)) <= f2c_dabs(tol) * (d__1 = d__[m], f2c_dabs(d__1))
           || tol < 0. && (d__3 = e[m - 1], f2c_dabs(d__3)) <= thresh)
        {
            e[m - 1] = 0.;
            goto L60;
        }
        if(tol >= 0.)
        {
            /* If relative accuracy desired, */
            /* apply convergence criterion forward */
            mu = (d__1 = d__[ll], f2c_dabs(d__1));
            smin = mu;
            i__1 = m - 1;
            for(lll = ll; lll <= i__1; ++lll)
            {
                if((d__1 = e[lll], f2c_dabs(d__1)) <= tol * mu)
                {
                    e[lll] = 0.;
                    goto L60;
                }
                mu = (d__2 = d__[lll + 1], f2c_dabs(d__2))
                     * (mu / (mu + (d__1 = e[lll], f2c_dabs(d__1))));
                smin = fla_min(smin, mu);
                /* L100: */
            }
        }
    }
    else
    {
        /* Run convergence test in backward direction */
        /* First apply standard test to top of matrix */
        if((d__2 = e[ll], f2c_dabs(d__2)) <= f2c_dabs(tol) * (d__1 = d__[ll], f2c_dabs(d__1))
           || tol < 0. && (d__3 = e[ll], f2c_dabs(d__3)) <= thresh)
        {
            e[ll] = 0.;
            goto L60;
        }
        if(tol >= 0.)
        {
            /* If relative accuracy desired, */
            /* apply convergence criterion backward */
            mu = (d__1 = d__[m], f2c_dabs(d__1));
            smin = mu;
            i__1 = ll;
            for(lll = m - 1; lll >= i__1; --lll)
            {
                if((d__1 = e[lll], f2c_dabs(d__1)) <= tol * mu)
                {
                    e[lll] = 0.;
                    goto L60;
                }
                mu = (d__2 = d__[lll], f2c_dabs(d__2))
                     * (mu / (mu + (d__1 = e[lll], f2c_dabs(d__1))));
                smin = fla_min(smin, mu);
                /* L110: */
            }
        }
    }
    oldll = ll;
    oldm = m;
    /* Compute shift. First, test if shifting would ruin relative */
    /* accuracy, and if so set the shift to zero. */
    /* Computing MAX */
    d__1 = eps;
    d__2 = tol * .01; // , expr subst
    if(tol >= 0. && *n * tol * (smin / smax) <= fla_max(d__1, d__2))
    {
        /* Use a zero shift to avoid loss of relative accuracy */
        shift = 0.;
    }
    else
    {
        /* Compute the shift from 2-by-2 block at end of matrix */
        if(idir == 1)
        {
            sll = (d__1 = d__[ll], f2c_dabs(d__1));
            dlas2_(&d__[m - 1], &e[m - 1], &d__[m], &shift, &r__);
        }
        else
        {
            sll = (d__1 = d__[m], f2c_dabs(d__1));
            dlas2_(&d__[ll], &e[ll], &d__[ll + 1], &shift, &r__);
        }
        /* Test if shift negligible, and if so set to zero */
        if(sll > 0.)
        {
            /* Computing 2nd power */
            d__1 = shift / sll;
            if(d__1 * d__1 < eps)
            {
                shift = 0.;
            }
        }
    }
    /* Increment iteration count */
    iter = iter + m - ll;
    /* If SHIFT = 0, do simplified QR iteration */
    if(shift == 0.)
    {
        if(idir == 1)
        {
            /* Chase bulge from top to bottom */
            /* Save cosines and sines for later singular vector updates */
            cs = 1.;
            oldcs = 1.;
            i__1 = m - 1;
            for(i__ = ll; i__ <= i__1; ++i__)
            {
                d__1 = d__[i__] * cs;
                dlartg_(&d__1, &e[i__], &cs, &sn, &r__);
                if(i__ > ll)
                {
                    e[i__ - 1] = oldsn * r__;
                }
                d__1 = oldcs * r__;
                d__2 = d__[i__ + 1] * sn;
                dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
                rwork[i__ - ll + 1] = cs;
                rwork[i__ - ll + 1 + nm1] = sn;
                rwork[i__ - ll + 1 + nm12] = oldcs;
                rwork[i__ - ll + 1 + nm13] = oldsn;
                /* L120: */
            }
            h__ = d__[m] * cs;
            d__[m] = h__ * oldcs;
            e[m - 1] = h__ * oldsn;
            /* Update singular vectors */
            if(*ncvt > 0)
            {
                i__1 = m - ll + 1;
                zlasr_("L", "V", "F", &i__1, ncvt, &rwork[1], &rwork[*n], &vt[ll + vt_dim1], ldvt);
            }
            if(*nru > 0)
            {
                i__1 = m - ll + 1;
                zlasr_("R", "V", "F", nru, &i__1, &rwork[nm12 + 1], &rwork[nm13 + 1],
                       &u[ll * u_dim1 + 1], ldu);
            }
            if(*ncc > 0)
            {
                i__1 = m - ll + 1;
                zlasr_("L", "V", "F", &i__1, ncc, &rwork[nm12 + 1], &rwork[nm13 + 1],
                       &c__[ll + c_dim1], ldc);
            }
            /* Test convergence */
            if((d__1 = e[m - 1], f2c_dabs(d__1)) <= thresh)
            {
                e[m - 1] = 0.;
            }
        }
        else
        {
            /* Chase bulge from bottom to top */
            /* Save cosines and sines for later singular vector updates */
            cs = 1.;
            oldcs = 1.;
            i__1 = ll + 1;
            for(i__ = m; i__ >= i__1; --i__)
            {
                d__1 = d__[i__] * cs;
                dlartg_(&d__1, &e[i__ - 1], &cs, &sn, &r__);
                if(i__ < m)
                {
                    e[i__] = oldsn * r__;
                }
                d__1 = oldcs * r__;
                d__2 = d__[i__ - 1] * sn;
                dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
                rwork[i__ - ll] = cs;
                rwork[i__ - ll + nm1] = -sn;
                rwork[i__ - ll + nm12] = oldcs;
                rwork[i__ - ll + nm13] = -oldsn;
                /* L130: */
            }
            h__ = d__[ll] * cs;
            d__[ll] = h__ * oldcs;
            e[ll] = h__ * oldsn;
            /* Update singular vectors */
            if(*ncvt > 0)
            {
                i__1 = m - ll + 1;
                zlasr_("L", "V", "B", &i__1, ncvt, &rwork[nm12 + 1], &rwork[nm13 + 1],
                       &vt[ll + vt_dim1], ldvt);
            }
            if(*nru > 0)
            {
                i__1 = m - ll + 1;
                zlasr_("R", "V", "B", nru, &i__1, &rwork[1], &rwork[*n], &u[ll * u_dim1 + 1], ldu);
            }
            if(*ncc > 0)
            {
                i__1 = m - ll + 1;
                zlasr_("L", "V", "B", &i__1, ncc, &rwork[1], &rwork[*n], &c__[ll + c_dim1], ldc);
            }
            /* Test convergence */
            if((d__1 = e[ll], f2c_dabs(d__1)) <= thresh)
            {
                e[ll] = 0.;
            }
        }
    }
    else
    {
        /* Use nonzero shift */
        if(idir == 1)
        {
            /* Chase bulge from top to bottom */
            /* Save cosines and sines for later singular vector updates */
            f = ((d__1 = d__[ll], f2c_dabs(d__1)) - shift)
                * (d_sign(&c_b49, &d__[ll]) + shift / d__[ll]);
            g = e[ll];
            i__1 = m - 1;
            for(i__ = ll; i__ <= i__1; ++i__)
            {
                dlartg_(&f, &g, &cosr, &sinr, &r__);
                if(i__ > ll)
                {
                    e[i__ - 1] = r__;
                }
                f = cosr * d__[i__] + sinr * e[i__];
                e[i__] = cosr * e[i__] - sinr * d__[i__];
                g = sinr * d__[i__ + 1];
                d__[i__ + 1] = cosr * d__[i__ + 1];
                dlartg_(&f, &g, &cosl, &sinl, &r__);
                d__[i__] = r__;
                f = cosl * e[i__] + sinl * d__[i__ + 1];
                d__[i__ + 1] = cosl * d__[i__ + 1] - sinl * e[i__];
                if(i__ < m - 1)
                {
                    g = sinl * e[i__ + 1];
                    e[i__ + 1] = cosl * e[i__ + 1];
                }
                rwork[i__ - ll + 1] = cosr;
                rwork[i__ - ll + 1 + nm1] = sinr;
                rwork[i__ - ll + 1 + nm12] = cosl;
                rwork[i__ - ll + 1 + nm13] = sinl;
                /* L140: */
            }
            e[m - 1] = f;
            /* Update singular vectors */
            if(*ncvt > 0)
            {
                i__1 = m - ll + 1;
                zlasr_("L", "V", "F", &i__1, ncvt, &rwork[1], &rwork[*n], &vt[ll + vt_dim1], ldvt);
            }
            if(*nru > 0)
            {
                i__1 = m - ll + 1;
                zlasr_("R", "V", "F", nru, &i__1, &rwork[nm12 + 1], &rwork[nm13 + 1],
                       &u[ll * u_dim1 + 1], ldu);
            }
            if(*ncc > 0)
            {
                i__1 = m - ll + 1;
                zlasr_("L", "V", "F", &i__1, ncc, &rwork[nm12 + 1], &rwork[nm13 + 1],
                       &c__[ll + c_dim1], ldc);
            }
            /* Test convergence */
            if((d__1 = e[m - 1], f2c_dabs(d__1)) <= thresh)
            {
                e[m - 1] = 0.;
            }
        }
        else
        {
            /* Chase bulge from bottom to top */
            /* Save cosines and sines for later singular vector updates */
            f = ((d__1 = d__[m], f2c_dabs(d__1)) - shift)
                * (d_sign(&c_b49, &d__[m]) + shift / d__[m]);
            g = e[m - 1];
            i__1 = ll + 1;
            for(i__ = m; i__ >= i__1; --i__)
            {
                dlartg_(&f, &g, &cosr, &sinr, &r__);
                if(i__ < m)
                {
                    e[i__] = r__;
                }
                f = cosr * d__[i__] + sinr * e[i__ - 1];
                e[i__ - 1] = cosr * e[i__ - 1] - sinr * d__[i__];
                g = sinr * d__[i__ - 1];
                d__[i__ - 1] = cosr * d__[i__ - 1];
                dlartg_(&f, &g, &cosl, &sinl, &r__);
                d__[i__] = r__;
                f = cosl * e[i__ - 1] + sinl * d__[i__ - 1];
                d__[i__ - 1] = cosl * d__[i__ - 1] - sinl * e[i__ - 1];
                if(i__ > ll + 1)
                {
                    g = sinl * e[i__ - 2];
                    e[i__ - 2] = cosl * e[i__ - 2];
                }
                rwork[i__ - ll] = cosr;
                rwork[i__ - ll + nm1] = -sinr;
                rwork[i__ - ll + nm12] = cosl;
                rwork[i__ - ll + nm13] = -sinl;
                /* L150: */
            }
            e[ll] = f;
            /* Test convergence */
            if((d__1 = e[ll], f2c_dabs(d__1)) <= thresh)
            {
                e[ll] = 0.;
            }
            /* Update singular vectors if desired */
            if(*ncvt > 0)
            {
                i__1 = m - ll + 1;
                zlasr_("L", "V", "B", &i__1, ncvt, &rwork[nm12 + 1], &rwork[nm13 + 1],
                       &vt[ll + vt_dim1], ldvt);
            }
            if(*nru > 0)
            {
                i__1 = m - ll + 1;
                zlasr_("R", "V", "B", nru, &i__1, &rwork[1], &rwork[*n], &u[ll * u_dim1 + 1], ldu);
            }
            if(*ncc > 0)
            {
                i__1 = m - ll + 1;
                zlasr_("L", "V", "B", &i__1, ncc, &rwork[1], &rwork[*n], &c__[ll + c_dim1], ldc);
            }
        }
    }
    /* QR iteration finished, go back and check convergence */
    goto L60;
/* All singular values converged, so make them positive */
L160:
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if(d__[i__] < 0.)
        {
            d__[i__] = -d__[i__];
            /* Change sign of singular vectors, if desired */
            if(*ncvt > 0)
            {
                zdscal_(ncvt, &c_b72, &vt[i__ + vt_dim1], ldvt);
            }
        }
        /* L170: */
    }
    /* Sort the singular values into decreasing order (insertion sort on */
    /* singular values, but only one transposition per singular vector) */
    i__1 = *n - 1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        /* Scan for smallest D(I) */
        isub = 1;
        smin = d__[1];
        i__2 = *n + 1 - i__;
        for(j = 2; j <= i__2; ++j)
        {
            if(d__[j] <= smin)
            {
                isub = j;
                smin = d__[j];
            }
            /* L180: */
        }
        if(isub != *n + 1 - i__)
        {
            /* Swap singular values and vectors */
            d__[isub] = d__[*n + 1 - i__];
            d__[*n + 1 - i__] = smin;
            if(*ncvt > 0)
            {
                zswap_(ncvt, &vt[isub + vt_dim1], ldvt, &vt[*n + 1 - i__ + vt_dim1], ldvt);
            }
            if(*nru > 0)
            {
                zswap_(nru, &u[isub * u_dim1 + 1], &c__1, &u[(*n + 1 - i__) * u_dim1 + 1], &c__1);
            }
            if(*ncc > 0)
            {
                zswap_(ncc, &c__[isub + c_dim1], ldc, &c__[*n + 1 - i__ + c_dim1], ldc);
            }
        }
        /* L190: */
    }
    goto L220;
/* Maximum number of iterations exceeded, failure to converge */
L200:
    *info = 0;
    i__1 = *n - 1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if(e[i__] != 0.)
        {
            ++(*info);
        }
        /* L210: */
    }
L220:
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZBDSQR */
}
/* zbdsqr_ */
