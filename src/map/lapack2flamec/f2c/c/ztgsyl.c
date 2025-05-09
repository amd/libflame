/* ../netlib/ztgsyl.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 = {0., 0.};
static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__5 = 5;
static integer c__1 = 1;
static doublecomplex c_b44 = {-1., 0.};
static doublecomplex c_b45 = {1., 0.};
/* > \brief \b ZTGSYL */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZTGSYL + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsyl.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsyl.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsyl.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, */
/* LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK, */
/* IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, */
/* $ LWORK, M, N */
/* DOUBLE PRECISION DIF, SCALE */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* COMPLEX*16 A( LDA, * ), B( LDB, * ), C( LDC, * ), */
/* $ D( LDD, * ), E( LDE, * ), F( LDF, * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTGSYL solves the generalized Sylvester equation: */
/* > */
/* > A * R - L * B = scale * C (1) */
/* > D * R - L * E = scale * F */
/* > */
/* > where R and L are unknown m-by-n matrices, (A, D), (B, E) and */
/* > (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n, */
/* > respectively, with complex entries. A, B, D and E are upper */
/* > triangular (i.e., (A,D) and (B,E) in generalized Schur form). */
/* > */
/* > The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 */
/* > is an output scaling factor chosen to avoid overflow. */
/* > */
/* > In matrix notation (1) is equivalent to solve Zx = scale*b, where Z */
/* > is defined as */
/* > */
/* > Z = [ kron(In, A) -kron(B**H, Im) ] (2) */
/* > [ kron(In, D) -kron(E**H, Im) ], */
/* > */
/* > Here Ix is the identity matrix of size x and X**H is the conjugate */
/* > transpose of X. Kron(X, Y) is the Kronecker product between the */
/* > matrices X and Y. */
/* > */
/* > If TRANS = 'C', y in the conjugate transposed system Z**H *y = scale*b */
/* > is solved for, which is equivalent to solve for R and L in */
/* > */
/* > A**H * R + D**H * L = scale * C (3) */
/* > R * B**H + L * E**H = scale * -F */
/* > */
/* > This case (TRANS = 'C') is used to compute an one-norm-based estimate */
/* > of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D) */
/* > and (B,E), using ZLACON. */
/* > */
/* > If IJOB >= 1, ZTGSYL computes a Frobenius norm-based estimate of */
/* > Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the */
/* > reciprocal of the smallest singular value of Z. */
/* > */
/* > This is a level-3 BLAS algorithm. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N': solve the generalized sylvester equation (1). */
/* > = 'C': solve the "conjugate transposed" system (3). */
/* > \endverbatim */
/* > */
/* > \param[in] IJOB */
/* > \verbatim */
/* > IJOB is INTEGER */
/* > Specifies what kind of functionality to be performed. */
/* > =0: solve (1) only. */
/* > =1: The functionality of 0 and 3. */
/* > =2: The functionality of 0 and 4. */
/* > =3: Only an estimate of Dif[(A,D), (B,E)] is computed. */
/* > (look ahead strategy is used). */
/* > =4: Only an estimate of Dif[(A,D), (B,E)] is computed. */
/* > (ZGECON on sub-systems is used). */
/* > Not referenced if TRANS = 'C'. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The order of the matrices A and D, and the row dimension of */
/* > the matrices C, F, R and L. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices B and E, and the column dimension */
/* > of the matrices C, F, R and L. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA, M) */
/* > The upper triangular matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB, N) */
/* > The upper triangular matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX*16 array, dimension (LDC, N) */
/* > On entry, C contains the right-hand-side of the first matrix */
/* > equation in (1) or (3). */
/* > On exit, if IJOB = 0, 1 or 2, C has been overwritten by */
/* > the solution R. If IJOB = 3 or 4 and TRANS = 'N', C holds R, */
/* > the solution achieved during the computation of the */
/* > Dif-estimate. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= fla_max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is COMPLEX*16 array, dimension (LDD, M) */
/* > The upper triangular matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] LDD */
/* > \verbatim */
/* > LDD is INTEGER */
/* > The leading dimension of the array D. LDD >= fla_max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is COMPLEX*16 array, dimension (LDE, N) */
/* > The upper triangular matrix E. */
/* > \endverbatim */
/* > */
/* > \param[in] LDE */
/* > \verbatim */
/* > LDE is INTEGER */
/* > The leading dimension of the array E. LDE >= fla_max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] F */
/* > \verbatim */
/* > F is COMPLEX*16 array, dimension (LDF, N) */
/* > On entry, F contains the right-hand-side of the second matrix */
/* > equation in (1) or (3). */
/* > On exit, if IJOB = 0, 1 or 2, F has been overwritten by */
/* > the solution L. If IJOB = 3 or 4 and TRANS = 'N', F holds L, */
/* > the solution achieved during the computation of the */
/* > Dif-estimate. */
/* > \endverbatim */
/* > */
/* > \param[in] LDF */
/* > \verbatim */
/* > LDF is INTEGER */
/* > The leading dimension of the array F. LDF >= fla_max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[out] DIF */
/* > \verbatim */
/* > DIF is DOUBLE PRECISION */
/* > On exit DIF is the reciprocal of a lower bound of the */
/* > reciprocal of the Dif-function, i.e. DIF is an upper bound of */
/* > Dif[(A,D), (B,E)] = sigma-min(Z), where Z as in (2). */
/* > IF IJOB = 0 or TRANS = 'C', DIF is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* > SCALE is DOUBLE PRECISION */
/* > On exit SCALE is the scaling factor in (1) or (3). */
/* > If 0 < SCALE < 1, C and F hold the solutions R and L, resp., */
/* > to a slightly perturbed system but the input matrices A, B, */
/* > D and E have not been changed. If SCALE = 0, R and L will */
/* > hold the solutions to the homogenious system with C = F = 0. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK > = 1. */
/* > If IJOB = 1 or 2 and TRANS = 'N', LWORK >= fla_max(1,2*M*N). */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (M+N+2) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > =0: successful exit */
/* > <0: If INFO = -i, the i-th argument had an illegal value. */
/* > >0: (A, D) and (B, E) have common or very close */
/* > eigenvalues. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complex16SYcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* > Umea University, S-901 87 Umea, Sweden. */
/* > \par References: */
/* ================ */
/* > */
/* > [1] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software */
/* > for Solving the Generalized Sylvester Equation and Estimating the */
/* > Separation between Regular Matrix Pairs, Report UMINF - 93.23, */
/* > Department of Computing Science, Umea University, S-901 87 Umea, */
/* > Sweden, December 1993, Revised April 1994, Also as LAPACK Working */
/* > Note 75. To appear in ACM Trans. on Math. Software, Vol 22, */
/* > No 1, 1996. */
/* > \n */
/* > [2] B. Kagstrom, A Perturbation Analysis of the Generalized Sylvester */
/* > Equation (AR - LB, DR - LE ) = (C, F), SIAM J. Matrix Anal. */
/* > Appl., 15(4):1045-1060, 1994. */
/* > \n */
/* > [3] B. Kagstrom and L. Westin, Generalized Schur Methods with */
/* > Condition Estimators for Solving the Generalized Sylvester */
/* > Equation, IEEE Transactions on Automatic Control, Vol. 34, No. 7, */
/* > July 1989, pp 745-751. */
/* > */
/* ===================================================================== */
/* Subroutine */
void ztgsyl_(char *trans, integer *ijob, integer *m, integer *n, doublecomplex *a, integer *lda,
             doublecomplex *b, integer *ldb, doublecomplex *c__, integer *ldc, doublecomplex *d__,
             integer *ldd, doublecomplex *e, integer *lde, doublecomplex *f, integer *ldf,
             doublereal *scale, doublereal *dif, doublecomplex *work, integer *lwork,
             integer *iwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("ztgsyl inputs: trans %c, ijob %" FLA_IS ", m %" FLA_IS ", n %" FLA_IS
                      ", lda %" FLA_IS ", ldb %" FLA_IS ", ldc %" FLA_IS ", ldd %" FLA_IS
                      ", lde %" FLA_IS ", ldf %" FLA_IS "",
                      *trans, *ijob, *m, *n, *lda, *ldb, *ldc, *ldd, *lde, *ldf);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, d_offset, e_dim1,
        e_offset, f_dim1, f_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, k, p, q, ie, je, mb, nb, is, js, pq;
    doublereal dsum;
    extern logical lsame_(char *, char *, integer, integer);
    integer ifunc, linfo;
    extern /* Subroutine */
        void
        zscal_(integer *, doublecomplex *, doublecomplex *, integer *),
        zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *, doublecomplex *,
               integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *);
    integer lwmin;
    doublereal scale2, dscale;
    extern /* Subroutine */
        void
        ztgsy2_(char *, integer *, integer *, integer *, doublecomplex *, integer *,
                doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *,
                doublecomplex *, integer *, doublecomplex *, integer *, doublereal *, doublereal *,
                doublereal *, integer *);
    doublereal scaloc;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer iround;
    logical notran;
    integer isolve;
    extern /* Subroutine */
        void
        zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *,
                integer *),
        zlaset_(char *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *,
                integer *);
    logical lquery;
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* Replaced various illegal calls to CCOPY by calls to CLASET. */
    /* Sven Hammarling, 1/5/02. */
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
    /* Decode and test input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    notran = lsame_(trans, "N", 1, 1);
    lquery = *lwork == -1;
    scale2 = 0.;
    if(!notran && !lsame_(trans, "C", 1, 1))
    {
        *info = -1;
    }
    else if(notran)
    {
        if(*ijob < 0 || *ijob > 4)
        {
            *info = -2;
        }
    }
    if(*info == 0)
    {
        if(*m <= 0)
        {
            *info = -3;
        }
        else if(*n <= 0)
        {
            *info = -4;
        }
        else if(*lda < fla_max(1, *m))
        {
            *info = -6;
        }
        else if(*ldb < fla_max(1, *n))
        {
            *info = -8;
        }
        else if(*ldc < fla_max(1, *m))
        {
            *info = -10;
        }
        else if(*ldd < fla_max(1, *m))
        {
            *info = -12;
        }
        else if(*lde < fla_max(1, *n))
        {
            *info = -14;
        }
        else if(*ldf < fla_max(1, *m))
        {
            *info = -16;
        }
    }
    if(*info == 0)
    {
        if(notran)
        {
            if(*ijob == 1 || *ijob == 2)
            {
                /* Computing MAX */
                i__1 = 1;
                i__2 = (*m << 1) * *n; // , expr subst
                lwmin = fla_max(i__1, i__2);
            }
            else
            {
                lwmin = 1;
            }
        }
        else
        {
            lwmin = 1;
        }
        work[1].r = (doublereal)lwmin;
        work[1].i = 0.; // , expr subst
        if(*lwork < lwmin && !lquery)
        {
            *info = -20;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZTGSYL", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*m == 0 || *n == 0)
    {
        *scale = 1.;
        if(notran)
        {
            if(*ijob != 0)
            {
                *dif = 0.;
            }
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Determine optimal block sizes MB and NB */
    mb = ilaenv_(&c__2, "ZTGSYL", trans, m, n, &c_n1, &c_n1);
    nb = ilaenv_(&c__5, "ZTGSYL", trans, m, n, &c_n1, &c_n1);
    isolve = 1;
    ifunc = 0;
    if(notran)
    {
        if(*ijob >= 3)
        {
            ifunc = *ijob - 2;
            zlaset_("F", m, n, &c_b1, &c_b1, &c__[c_offset], ldc);
            zlaset_("F", m, n, &c_b1, &c_b1, &f[f_offset], ldf);
        }
        else if(*ijob >= 1 && notran)
        {
            isolve = 2;
        }
    }
    if(mb <= 1 && nb <= 1 || mb >= *m && nb >= *n)
    {
        /* Use unblocked Level 2 solver */
        i__1 = isolve;
        for(iround = 1; iround <= i__1; ++iround)
        {
            *scale = 1.;
            dscale = 0.;
            dsum = 1.;
            pq = *m * *n;
            ztgsy2_(trans, &ifunc, m, n, &a[a_offset], lda, &b[b_offset], ldb, &c__[c_offset], ldc,
                    &d__[d_offset], ldd, &e[e_offset], lde, &f[f_offset], ldf, scale, &dsum,
                    &dscale, info);
            if(dscale != 0.)
            {
                if(*ijob == 1 || *ijob == 3)
                {
                    *dif = sqrt((doublereal)((*m << 1) * *n)) / (dscale * sqrt(dsum));
                }
                else
                {
                    *dif = sqrt((doublereal)pq) / (dscale * sqrt(dsum));
                }
            }
            if(isolve == 2 && iround == 1)
            {
                if(notran)
                {
                    ifunc = *ijob;
                }
                scale2 = *scale;
                zlacpy_("F", m, n, &c__[c_offset], ldc, &work[1], m);
                zlacpy_("F", m, n, &f[f_offset], ldf, &work[*m * *n + 1], m);
                zlaset_("F", m, n, &c_b1, &c_b1, &c__[c_offset], ldc);
                zlaset_("F", m, n, &c_b1, &c_b1, &f[f_offset], ldf);
            }
            else if(isolve == 2 && iround == 2)
            {
                zlacpy_("F", m, n, &work[1], m, &c__[c_offset], ldc);
                zlacpy_("F", m, n, &work[*m * *n + 1], m, &f[f_offset], ldf);
                *scale = scale2;
            }
            /* L30: */
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Determine block structure of A */
    p = 0;
    i__ = 1;
L40:
    if(i__ > *m)
    {
        goto L50;
    }
    ++p;
    iwork[p] = i__;
    i__ += mb;
    if(i__ >= *m)
    {
        goto L50;
    }
    goto L40;
L50:
    iwork[p + 1] = *m + 1;
    if(iwork[p] == iwork[p + 1])
    {
        --p;
    }
    /* Determine block structure of B */
    q = p + 1;
    j = 1;
L60:
    if(j > *n)
    {
        goto L70;
    }
    ++q;
    iwork[q] = j;
    j += nb;
    if(j >= *n)
    {
        goto L70;
    }
    goto L60;
L70:
    iwork[q + 1] = *n + 1;
    if(iwork[q] == iwork[q + 1])
    {
        --q;
    }
    if(notran)
    {
        i__1 = isolve;
        for(iround = 1; iround <= i__1; ++iround)
        {
            /* Solve (I, J) - subsystem */
            /* A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J) */
            /* D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J) */
            /* for I = P, P - 1, ..., 1;
            J = 1, 2, ..., Q */
            pq = 0;
            *scale = 1.;
            dscale = 0.;
            dsum = 1.;
            i__2 = q;
            for(j = p + 2; j <= i__2; ++j)
            {
                js = iwork[j];
                je = iwork[j + 1] - 1;
                nb = je - js + 1;
                for(i__ = p; i__ >= 1; --i__)
                {
                    is = iwork[i__];
                    ie = iwork[i__ + 1] - 1;
                    mb = ie - is + 1;
                    ztgsy2_(trans, &ifunc, &mb, &nb, &a[is + is * a_dim1], lda,
                            &b[js + js * b_dim1], ldb, &c__[is + js * c_dim1], ldc,
                            &d__[is + is * d_dim1], ldd, &e[js + js * e_dim1], lde,
                            &f[is + js * f_dim1], ldf, &scaloc, &dsum, &dscale, &linfo);
                    if(linfo > 0)
                    {
                        *info = linfo;
                    }
                    pq += mb * nb;
                    if(scaloc != 1.)
                    {
                        i__3 = js - 1;
                        for(k = 1; k <= i__3; ++k)
                        {
                            z__1.r = scaloc;
                            z__1.i = 0.; // , expr subst
                            zscal_(m, &z__1, &c__[k * c_dim1 + 1], &c__1);
                            z__1.r = scaloc;
                            z__1.i = 0.; // , expr subst
                            zscal_(m, &z__1, &f[k * f_dim1 + 1], &c__1);
                            /* L80: */
                        }
                        i__3 = je;
                        for(k = js; k <= i__3; ++k)
                        {
                            i__4 = is - 1;
                            z__1.r = scaloc;
                            z__1.i = 0.; // , expr subst
                            zscal_(&i__4, &z__1, &c__[k * c_dim1 + 1], &c__1);
                            i__4 = is - 1;
                            z__1.r = scaloc;
                            z__1.i = 0.; // , expr subst
                            zscal_(&i__4, &z__1, &f[k * f_dim1 + 1], &c__1);
                            /* L90: */
                        }
                        i__3 = je;
                        for(k = js; k <= i__3; ++k)
                        {
                            i__4 = *m - ie;
                            z__1.r = scaloc;
                            z__1.i = 0.; // , expr subst
                            zscal_(&i__4, &z__1, &c__[ie + 1 + k * c_dim1], &c__1);
                            i__4 = *m - ie;
                            z__1.r = scaloc;
                            z__1.i = 0.; // , expr subst
                            zscal_(&i__4, &z__1, &f[ie + 1 + k * f_dim1], &c__1);
                            /* L100: */
                        }
                        i__3 = *n;
                        for(k = je + 1; k <= i__3; ++k)
                        {
                            z__1.r = scaloc;
                            z__1.i = 0.; // , expr subst
                            zscal_(m, &z__1, &c__[k * c_dim1 + 1], &c__1);
                            z__1.r = scaloc;
                            z__1.i = 0.; // , expr subst
                            zscal_(m, &z__1, &f[k * f_dim1 + 1], &c__1);
                            /* L110: */
                        }
                        *scale *= scaloc;
                    }
                    /* Substitute R(I,J) and L(I,J) into remaining equation. */
                    if(i__ > 1)
                    {
                        i__3 = is - 1;
                        zgemm_("N", "N", &i__3, &nb, &mb, &c_b44, &a[is * a_dim1 + 1], lda,
                               &c__[is + js * c_dim1], ldc, &c_b45, &c__[js * c_dim1 + 1], ldc);
                        i__3 = is - 1;
                        zgemm_("N", "N", &i__3, &nb, &mb, &c_b44, &d__[is * d_dim1 + 1], ldd,
                               &c__[is + js * c_dim1], ldc, &c_b45, &f[js * f_dim1 + 1], ldf);
                    }
                    if(j < q)
                    {
                        i__3 = *n - je;
                        zgemm_("N", "N", &mb, &i__3, &nb, &c_b45, &f[is + js * f_dim1], ldf,
                               &b[js + (je + 1) * b_dim1], ldb, &c_b45,
                               &c__[is + (je + 1) * c_dim1], ldc);
                        i__3 = *n - je;
                        zgemm_("N", "N", &mb, &i__3, &nb, &c_b45, &f[is + js * f_dim1], ldf,
                               &e[js + (je + 1) * e_dim1], lde, &c_b45, &f[is + (je + 1) * f_dim1],
                               ldf);
                    }
                    /* L120: */
                }
                /* L130: */
            }
            if(dscale != 0.)
            {
                if(*ijob == 1 || *ijob == 3)
                {
                    *dif = sqrt((doublereal)((*m << 1) * *n)) / (dscale * sqrt(dsum));
                }
                else
                {
                    *dif = sqrt((doublereal)pq) / (dscale * sqrt(dsum));
                }
            }
            if(isolve == 2 && iround == 1)
            {
                if(notran)
                {
                    ifunc = *ijob;
                }
                scale2 = *scale;
                zlacpy_("F", m, n, &c__[c_offset], ldc, &work[1], m);
                zlacpy_("F", m, n, &f[f_offset], ldf, &work[*m * *n + 1], m);
                zlaset_("F", m, n, &c_b1, &c_b1, &c__[c_offset], ldc);
                zlaset_("F", m, n, &c_b1, &c_b1, &f[f_offset], ldf);
            }
            else if(isolve == 2 && iround == 2)
            {
                zlacpy_("F", m, n, &work[1], m, &c__[c_offset], ldc);
                zlacpy_("F", m, n, &work[*m * *n + 1], m, &f[f_offset], ldf);
                *scale = scale2;
            }
            /* L150: */
        }
    }
    else
    {
        /* Solve transposed (I, J)-subsystem */
        /* A(I, I)**H * R(I, J) + D(I, I)**H * L(I, J) = C(I, J) */
        /* R(I, J) * B(J, J) + L(I, J) * E(J, J) = -F(I, J) */
        /* for I = 1,2,..., P;
        J = Q, Q-1,..., 1 */
        *scale = 1.;
        i__1 = p;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            is = iwork[i__];
            ie = iwork[i__ + 1] - 1;
            mb = ie - is + 1;
            i__2 = p + 2;
            for(j = q; j >= i__2; --j)
            {
                js = iwork[j];
                je = iwork[j + 1] - 1;
                nb = je - js + 1;
                ztgsy2_(trans, &ifunc, &mb, &nb, &a[is + is * a_dim1], lda, &b[js + js * b_dim1],
                        ldb, &c__[is + js * c_dim1], ldc, &d__[is + is * d_dim1], ldd,
                        &e[js + js * e_dim1], lde, &f[is + js * f_dim1], ldf, &scaloc, &dsum,
                        &dscale, &linfo);
                if(linfo > 0)
                {
                    *info = linfo;
                }
                if(scaloc != 1.)
                {
                    i__3 = js - 1;
                    for(k = 1; k <= i__3; ++k)
                    {
                        z__1.r = scaloc;
                        z__1.i = 0.; // , expr subst
                        zscal_(m, &z__1, &c__[k * c_dim1 + 1], &c__1);
                        z__1.r = scaloc;
                        z__1.i = 0.; // , expr subst
                        zscal_(m, &z__1, &f[k * f_dim1 + 1], &c__1);
                        /* L160: */
                    }
                    i__3 = je;
                    for(k = js; k <= i__3; ++k)
                    {
                        i__4 = is - 1;
                        z__1.r = scaloc;
                        z__1.i = 0.; // , expr subst
                        zscal_(&i__4, &z__1, &c__[k * c_dim1 + 1], &c__1);
                        i__4 = is - 1;
                        z__1.r = scaloc;
                        z__1.i = 0.; // , expr subst
                        zscal_(&i__4, &z__1, &f[k * f_dim1 + 1], &c__1);
                        /* L170: */
                    }
                    i__3 = je;
                    for(k = js; k <= i__3; ++k)
                    {
                        i__4 = *m - ie;
                        z__1.r = scaloc;
                        z__1.i = 0.; // , expr subst
                        zscal_(&i__4, &z__1, &c__[ie + 1 + k * c_dim1], &c__1);
                        i__4 = *m - ie;
                        z__1.r = scaloc;
                        z__1.i = 0.; // , expr subst
                        zscal_(&i__4, &z__1, &f[ie + 1 + k * f_dim1], &c__1);
                        /* L180: */
                    }
                    i__3 = *n;
                    for(k = je + 1; k <= i__3; ++k)
                    {
                        z__1.r = scaloc;
                        z__1.i = 0.; // , expr subst
                        zscal_(m, &z__1, &c__[k * c_dim1 + 1], &c__1);
                        z__1.r = scaloc;
                        z__1.i = 0.; // , expr subst
                        zscal_(m, &z__1, &f[k * f_dim1 + 1], &c__1);
                        /* L190: */
                    }
                    *scale *= scaloc;
                }
                /* Substitute R(I,J) and L(I,J) into remaining equation. */
                if(j > p + 2)
                {
                    i__3 = js - 1;
                    zgemm_("N", "C", &mb, &i__3, &nb, &c_b45, &c__[is + js * c_dim1], ldc,
                           &b[js * b_dim1 + 1], ldb, &c_b45, &f[is + f_dim1], ldf);
                    i__3 = js - 1;
                    zgemm_("N", "C", &mb, &i__3, &nb, &c_b45, &f[is + js * f_dim1], ldf,
                           &e[js * e_dim1 + 1], lde, &c_b45, &f[is + f_dim1], ldf);
                }
                if(i__ < p)
                {
                    i__3 = *m - ie;
                    zgemm_("C", "N", &i__3, &nb, &mb, &c_b44, &a[is + (ie + 1) * a_dim1], lda,
                           &c__[is + js * c_dim1], ldc, &c_b45, &c__[ie + 1 + js * c_dim1], ldc);
                    i__3 = *m - ie;
                    zgemm_("C", "N", &i__3, &nb, &mb, &c_b44, &d__[is + (ie + 1) * d_dim1], ldd,
                           &f[is + js * f_dim1], ldf, &c_b45, &c__[ie + 1 + js * c_dim1], ldc);
                }
                /* L200: */
            }
            /* L210: */
        }
    }
    work[1].r = (doublereal)lwmin;
    work[1].i = 0.; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZTGSYL */
}
/* ztgsyl_ */
