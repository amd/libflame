/* ../netlib/zstein.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b ZSTEIN */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZSTEIN + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zstein.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zstein.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zstein.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, */
/* IWORK, IFAIL, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDZ, M, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IBLOCK( * ), IFAIL( * ), ISPLIT( * ), */
/* $ IWORK( * ) */
/* DOUBLE PRECISION D( * ), E( * ), W( * ), WORK( * ) */
/* COMPLEX*16 Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSTEIN computes the eigenvectors of a real symmetric tridiagonal */
/* > matrix T corresponding to specified eigenvalues, using inverse */
/* > iteration. */
/* > */
/* > The maximum number of iterations allowed for each eigenvector is */
/* > specified by an internal parameter MAXITS (currently set to 5). */
/* > */
/* > Although the eigenvectors are real, they are stored in a complex */
/* > array, which may be passed to ZUNMTR or ZUPMTR for back */
/* > transformation to the eigenvectors of a complex Hermitian matrix */
/* > which was reduced to tridiagonal form. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (N) */
/* > The n diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is DOUBLE PRECISION array, dimension (N-1) */
/* > The (n-1) subdiagonal elements of the tridiagonal matrix */
/* > T, stored in elements 1 to N-1. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of eigenvectors to be found. 0 <= M <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* > W is DOUBLE PRECISION array, dimension (N) */
/* > The first M elements of W contain the eigenvalues for */
/* > which eigenvectors are to be computed. The eigenvalues */
/* > should be grouped by split-off block and ordered from */
/* > smallest to largest within the block. ( The output array */
/* > W from DSTEBZ with ORDER = 'B' is expected here. ) */
/* > \endverbatim */
/* > */
/* > \param[in] IBLOCK */
/* > \verbatim */
/* > IBLOCK is INTEGER array, dimension (N) */
/* > The submatrix indices associated with the corresponding */
/* > eigenvalues in W;
IBLOCK(i)=1 if eigenvalue W(i) belongs to */
/* > the first submatrix from the top, =2 if W(i) belongs to */
/* > the second submatrix, etc. ( The output array IBLOCK */
/* > from DSTEBZ is expected here. ) */
/* > \endverbatim */
/* > */
/* > \param[in] ISPLIT */
/* > \verbatim */
/* > ISPLIT is INTEGER array, dimension (N) */
/* > The splitting points, at which T breaks up into submatrices. */
/* > The first submatrix consists of rows/columns 1 to */
/* > ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1 */
/* > through ISPLIT( 2 ), etc. */
/* > ( The output array ISPLIT from DSTEBZ is expected here. ) */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is COMPLEX*16 array, dimension (LDZ, M) */
/* > The computed eigenvectors. The eigenvector associated */
/* > with the eigenvalue W(i) is stored in the i-th column of */
/* > Z. Any vector which fails to converge is set to its current */
/* > iterate after MAXITS iterations. */
/* > The imaginary parts of the eigenvectors are set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (5*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] IFAIL */
/* > \verbatim */
/* > IFAIL is INTEGER array, dimension (M) */
/* > On normal exit, all elements of IFAIL are zero. */
/* > If one or more eigenvectors fail to converge after */
/* > MAXITS iterations, then their indices are stored in */
/* > array IFAIL. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, then i eigenvectors failed to converge */
/* > in MAXITS iterations. Their indices are stored in */
/* > array IFAIL. */
/* > \endverbatim */
/* > \par Internal Parameters: */
/* ========================= */
/* > */
/* > \verbatim */
/* > MAXITS INTEGER, default = 5 */
/* > The maximum number of iterations performed. */
/* > */
/* > EXTRA INTEGER, default = 2 */
/* > The number of iterations performed after norm growth */
/* > criterion is satisfied, should be at least 1. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup complex16OTHERcomputational */
/* ===================================================================== */
/* Subroutine */
void zstein_(integer *n, doublereal *d__, doublereal *e, integer *m, doublereal *w, integer *iblock,
             integer *isplit, doublecomplex *z__, integer *ldz, doublereal *work, integer *iwork,
             integer *ifail, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zstein inputs: n %" FLA_IS ", m %" FLA_IS ", ldz %" FLA_IS "", *n, *m, *ldz);
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5;
    doublecomplex z__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, b1, j1, bn, jr;
    doublereal xj, scl, eps, sep, nrm, tol;
    integer its;
    doublereal xjm, ztr, eps1;
    integer jblk, nblk, jmax;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */
        void
        dscal_(integer *, doublereal *, doublereal *, integer *);
    integer iseed[4], gpind, iinfo;
    extern /* Subroutine */
        void
        dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal ortol;
    integer indrv1, indrv2, indrv3, indrv4, indrv5;
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
        void
        dlagtf_(integer *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, integer *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern /* Subroutine */
        void
        dlagts_(integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                integer *, doublereal *, doublereal *, integer *);
    integer nrmchk;
    extern /* Subroutine */
        void
        dlarnv_(integer *, integer *, integer *, doublereal *);
    integer blksiz;
    doublereal onenrm, dtpcrt, pertol;
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
    --d__;
    --e;
    --w;
    --iblock;
    --isplit;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;
    --ifail;
    /* Function Body */
    *info = 0;
    i__1 = *m;
    dtpcrt = 0.;
    onenrm = 0.;
    ortol = 0.;
    gpind = 0;
    xjm = 0.;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        ifail[i__] = 0;
        /* L10: */
    }
    if(*n < 0)
    {
        *info = -1;
    }
    else if(*m < 0 || *m > *n)
    {
        *info = -4;
    }
    else if(*ldz < fla_max(1, *n))
    {
        *info = -9;
    }
    else
    {
        i__1 = *m;
        for(j = 2; j <= i__1; ++j)
        {
            if(iblock[j] < iblock[j - 1])
            {
                *info = -6;
                goto L30;
            }
            if(iblock[j] == iblock[j - 1] && w[j] < w[j - 1])
            {
                *info = -5;
                goto L30;
            }
            /* L20: */
        }
    L30:;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZSTEIN", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0 || *m == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(*n == 1)
    {
        i__1 = z_dim1 + 1;
        z__[i__1].r = 1.;
        z__[i__1].i = 0.; // , expr subst
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Get machine constants. */
    eps = dlamch_("Precision");
    /* Initialize seed for random number generator DLARNV. */
    for(i__ = 1; i__ <= 4; ++i__)
    {
        iseed[i__ - 1] = 1;
        /* L40: */
    }
    /* Initialize pointers. */
    indrv1 = 0;
    indrv2 = indrv1 + *n;
    indrv3 = indrv2 + *n;
    indrv4 = indrv3 + *n;
    indrv5 = indrv4 + *n;
    /* Compute eigenvectors of matrix blocks. */
    j1 = 1;
    i__1 = iblock[*m];
    for(nblk = 1; nblk <= i__1; ++nblk)
    {
        /* Find starting and ending indices of block nblk. */
        if(nblk == 1)
        {
            b1 = 1;
        }
        else
        {
            b1 = isplit[nblk - 1] + 1;
        }
        bn = isplit[nblk];
        blksiz = bn - b1 + 1;
        if(blksiz == 1)
        {
            goto L60;
        }
        gpind = j1;
        /* Compute reorthogonalization criterion and stopping criterion. */
        onenrm = (d__1 = d__[b1], f2c_abs(d__1)) + (d__2 = e[b1], f2c_abs(d__2));
        /* Computing MAX */
        d__3 = onenrm;
        d__4 = (d__1 = d__[bn], f2c_abs(d__1)) + (d__2 = e[bn - 1], f2c_abs(d__2)); // , expr subst
        onenrm = fla_max(d__3, d__4);
        i__2 = bn - 1;
        for(i__ = b1 + 1; i__ <= i__2; ++i__)
        {
            /* Computing MAX */
            d__4 = onenrm;
            d__5 = (d__1 = d__[i__], f2c_abs(d__1)) + (d__2 = e[i__ - 1], f2c_abs(d__2))
                   + (d__3 = e[i__], f2c_abs(d__3)); // , expr subst
            onenrm = fla_max(d__4, d__5);
            /* L50: */
        }
        ortol = onenrm * .001;
        dtpcrt = sqrt(.1 / blksiz);
        /* Loop through eigenvalues of block nblk. */
    L60:
        jblk = 0;
        i__2 = *m;
        for(j = j1; j <= i__2; ++j)
        {
            if(iblock[j] != nblk)
            {
                j1 = j;
                goto L180;
            }
            ++jblk;
            xj = w[j];
            /* Skip all the work if the block size is one. */
            if(blksiz == 1)
            {
                work[indrv1 + 1] = 1.;
                goto L140;
            }
            /* If eigenvalues j and j-1 are too close, add a relatively */
            /* small perturbation. */
            if(jblk > 1)
            {
                eps1 = (d__1 = eps * xj, f2c_abs(d__1));
                pertol = eps1 * 10.;
                sep = xj - xjm;
                if(sep < pertol)
                {
                    xj = xjm + pertol;
                }
            }
            its = 0;
            nrmchk = 0;
            /* Get random starting vector. */
            dlarnv_(&c__2, iseed, &blksiz, &work[indrv1 + 1]);
            /* Copy the matrix T so it won't be destroyed in factorization. */
            dcopy_(&blksiz, &d__[b1], &c__1, &work[indrv4 + 1], &c__1);
            i__3 = blksiz - 1;
            dcopy_(&i__3, &e[b1], &c__1, &work[indrv2 + 2], &c__1);
            i__3 = blksiz - 1;
            dcopy_(&i__3, &e[b1], &c__1, &work[indrv3 + 1], &c__1);
            /* Compute LU factors with partial pivoting ( PT = LU ) */
            tol = 0.;
            dlagtf_(&blksiz, &work[indrv4 + 1], &xj, &work[indrv2 + 2], &work[indrv3 + 1], &tol,
                    &work[indrv5 + 1], &iwork[1], &iinfo);
            /* Update iteration count. */
        L70:
            ++its;
            if(its > 5)
            {
                goto L120;
            }
            /* Normalize and scale the righthand side vector Pb. */
            jmax = idamax_(&blksiz, &work[indrv1 + 1], &c__1);
            /* Computing MAX */
            d__3 = eps;
            d__4 = (d__1 = work[indrv4 + blksiz], f2c_abs(d__1)); // , expr subst
            scl = blksiz * onenrm * fla_max(d__3, d__4)
                  / (d__2 = work[indrv1 + jmax], f2c_abs(d__2));
            dscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);
            /* Solve the system LU = Pb. */
            dlagts_(&c_n1, &blksiz, &work[indrv4 + 1], &work[indrv2 + 2], &work[indrv3 + 1],
                    &work[indrv5 + 1], &iwork[1], &work[indrv1 + 1], &tol, &iinfo);
            /* Reorthogonalize by modified Gram-Schmidt if eigenvalues are */
            /* close enough. */
            if(jblk == 1)
            {
                goto L110;
            }
            if((d__1 = xj - xjm, f2c_abs(d__1)) > ortol)
            {
                gpind = j;
            }
            if(gpind != j)
            {
                i__3 = j - 1;
                for(i__ = gpind; i__ <= i__3; ++i__)
                {
                    ztr = 0.;
                    i__4 = blksiz;
                    for(jr = 1; jr <= i__4; ++jr)
                    {
                        i__5 = b1 - 1 + jr + i__ * z_dim1;
                        ztr += work[indrv1 + jr] * z__[i__5].r;
                        /* L80: */
                    }
                    i__4 = blksiz;
                    for(jr = 1; jr <= i__4; ++jr)
                    {
                        i__5 = b1 - 1 + jr + i__ * z_dim1;
                        work[indrv1 + jr] -= ztr * z__[i__5].r;
                        /* L90: */
                    }
                    /* L100: */
                }
            }
            /* Check the infinity norm of the iterate. */
        L110:
            jmax = idamax_(&blksiz, &work[indrv1 + 1], &c__1);
            nrm = (d__1 = work[indrv1 + jmax], f2c_abs(d__1));
            /* Continue for additional iterations after norm reaches */
            /* stopping criterion. */
            if(nrm < dtpcrt)
            {
                goto L70;
            }
            ++nrmchk;
            if(nrmchk < 3)
            {
                goto L70;
            }
            goto L130;
            /* If stopping criterion was not satisfied, update info and */
            /* store eigenvector number in array ifail. */
        L120:
            ++(*info);
            ifail[*info] = j;
            /* Accept iterate as jth eigenvector. */
        L130:
            scl = 1. / dnrm2_(&blksiz, &work[indrv1 + 1], &c__1);
            jmax = idamax_(&blksiz, &work[indrv1 + 1], &c__1);
            if(work[indrv1 + jmax] < 0.)
            {
                scl = -scl;
            }
            dscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);
        L140:
            i__3 = *n;
            for(i__ = 1; i__ <= i__3; ++i__)
            {
                i__4 = i__ + j * z_dim1;
                z__[i__4].r = 0.;
                z__[i__4].i = 0.; // , expr subst
                /* L150: */
            }
            i__3 = blksiz;
            for(i__ = 1; i__ <= i__3; ++i__)
            {
                i__4 = b1 + i__ - 1 + j * z_dim1;
                i__5 = indrv1 + i__;
                z__1.r = work[i__5];
                z__1.i = 0.; // , expr subst
                z__[i__4].r = z__1.r;
                z__[i__4].i = z__1.i; // , expr subst
                /* L160: */
            }
            /* Save the shift to check eigenvalue spacing at next */
            /* iteration. */
            xjm = xj;
            /* L170: */
        }
    L180:;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZSTEIN */
}
/* zstein_ */
