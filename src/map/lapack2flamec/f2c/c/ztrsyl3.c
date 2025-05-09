/* ztrsyl3.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 = {1., 0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b18 = 2.;
static doublereal c_b106 = 1.;
/* > \brief \b ZTRSYL3 */
/* Definition: */
/* =========== */
/* > \par Purpose */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTRSYL3 solves the complex Sylvester matrix equation: */
/* > */
/* > op(A)*X + X*op(B) = scale*C or */
/* > op(A)*X - X*op(B) = scale*C, */
/* > */
/* > where op(A) = A or A**H, and A and B are both upper triangular. A is */
/* > M-by-M and B is N-by-N;
the right hand side C and the solution X are */
/* > M-by-N;
and scale is an output scale factor, set <= 1 to avoid */
/* > overflow in X. */
/* > */
/* > This is the block version of the algorithm. */
/* > \endverbatim */
/* Arguments */
/* ========= */
/* > \param[in] TRANA */
/* > \verbatim */
/* > TRANA is CHARACTER*1 */
/* > Specifies the option op(A): */
/* > = 'N': op(A) = A (No transpose) */
/* > = 'C': op(A) = A**H (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] TRANB */
/* > \verbatim */
/* > TRANB is CHARACTER*1 */
/* > Specifies the option op(B): */
/* > = 'N': op(B) = B (No transpose) */
/* > = 'C': op(B) = B**H (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] ISGN */
/* > \verbatim */
/* > ISGN is INTEGER */
/* > Specifies the sign in the equation: */
/* > = +1: solve op(A)*X + X*op(B) = scale*C */
/* > = -1: solve op(A)*X - X*op(B) = scale*C */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The order of the matrix A, and the number of rows in the */
/* > matrices X and C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix B, and the number of columns in the */
/* > matrices X and C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,M) */
/* > The upper triangular matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,N) */
/* > The upper triangular matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX*16 array, dimension (LDC,N) */
/* > On entry, the M-by-N right hand side matrix C. */
/* > On exit, C is overwritten by the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= fla_max(1,M) */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* > SCALE is DOUBLE PRECISION */
/* > The scale factor, scale, set <= 1 to avoid overflow in X. */
/* > \endverbatim */
/* > */
/* > \param[out] SWORK */
/* > \verbatim */
/* > SWORK is DOUBLE PRECISION array, dimension (fla_max(2, ROWS), */
/* > fla_max(1,COLS)). */
/* > On exit, if INFO = 0, SWORK(1) returns the optimal value ROWS */
/* > and SWORK(2) returns the optimal COLS. */
/* > \endverbatim */
/* > */
/* > \param[in] LDSWORK */
/* > \verbatim */
/* > LDSWORK is INTEGER */
/* > LDSWORK >= fla_max(2,ROWS), where ROWS = ((M + NB - 1) / NB + 1) */
/* > and NB is the optimal block size. */
/* > */
/* > If LDSWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal dimensions of the SWORK matrix, */
/* > returns these values as the first and second entry of the SWORK */
/* > matrix, and no error message related LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > = 1: A and B have common or very close eigenvalues;
perturbed */
/* > values were used to solve the equation (but the matrices */
/* > A and B are unchanged). */
/* > \endverbatim */
/* > \ingroup complex16SYcomputational */
/* ===================================================================== */
/* References: */
/* E. S. Quintana-Orti and R. A. Van De Geijn (2003). Formal derivation of */
/* algorithms: The triangular Sylvester equation, ACM Transactions */
/* on Mathematical Software (TOMS), volume 29, pages 218--243. */
/* A. Schwarz and C. C. Kjelgaard Mikkelsen (2020). Robust Task-Parallel */
/* Solution of the Triangular Sylvester Equation. Lecture Notes in */
/* Computer Science, vol 12043, pages 82--92, Springer. */
/* Contributor: */
/* Angelika Schwarz, Umea University, Sweden. */
/* ===================================================================== */
/* Subroutine */
void ztrsyl3_(char *trana, char *tranb, integer *isgn, integer *m, integer *n, doublecomplex *a,
              integer *lda, doublecomplex *b, integer *ldb, doublecomplex *c__, integer *ldc,
              doublereal *scale, doublereal *swork, integer *ldswork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("ztrsyl3 inputs: trana %c, tranb %c, isgn %" FLA_IS ", m %" FLA_IS
                      ", n %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS ", ldc %" FLA_IS "",
                      *trana, *tranb, *isgn, *m, *n, *lda, *ldb, *ldc);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, swork_dim1, swork_offset, i__1,
        i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;
    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), d_imag(doublecomplex *);
    /* Local variables */
    integer i__, j, k, l, i1, i2, j1, j2, k1, k2, l1, l2, nb, jj, ll, nba, nbb;
    doublereal buf, sgn, scal;
    doublecomplex csgn;
    doublereal anrm, bnrm, cnrm;
    integer awrk, bwrk;
    int temp;
    doublereal *wnrm, xnrm;
    extern logical lsame_(char *, char *, integer, integer);
    integer iinfo;
    extern /* Subroutine */
        void
        zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *, doublecomplex *,
               integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *);
    extern doublereal dlamch_(char *);
    doublereal scaloc, scamin;
    extern doublereal dlarmm_(doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, integer *,
                              doublereal *);
    doublereal bignum;
    extern /* Subroutine */
        void
        zdscal_(integer *, doublereal *, doublecomplex *, integer *),
        zlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublecomplex *, integer *, integer *);
    logical notrna, notrnb;
    doublereal smlnum;
    logical lquery;
    extern /* Subroutine */
        void
        ztrsyl_(char *, char *, integer *, integer *, integer *, doublecomplex *, integer *,
                doublecomplex *, integer *, doublecomplex *, integer *, doublereal *, integer *);
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Decode and Test input parameters */
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
    swork_dim1 = *ldswork;
    swork_offset = 1 + swork_dim1;
    swork -= swork_offset;
    wnrm = NULL;
    wnrm = (doublereal *)malloc(fla_max(*m, *n) * sizeof(doublereal));
    if(wnrm == NULL)
    {
        *info = -4;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Function Body */
    notrna = lsame_(trana, "N", 1, 1);
    notrnb = lsame_(tranb, "N", 1, 1);
    /* Use the same block size for all matrices. */
    /* Computing fla_max */
    i__1 = 8;
    i__2 = ilaenv_(&c__1, "ZTRSYL", "", m, n, &c_n1, &c_n1); // , expr subst
    nb = fla_max(i__1, i__2);
    /* Compute number of blocks in A and B */
    /* Computing fla_max */
    i__1 = 1;
    i__2 = (*m + nb - 1) / nb; // , expr subst
    nba = fla_max(i__1, i__2);
    /* Computing fla_max */
    i__1 = 1;
    i__2 = (*n + nb - 1) / nb; // , expr subst
    nbb = fla_max(i__1, i__2);
    /* Compute workspace */
    *info = 0;
    lquery = *ldswork == -1;
    if(lquery)
    {
        *ldswork = 2;
        swork[swork_dim1 + 1] = (doublereal)fla_max(nba, nbb);
        swork[swork_dim1 + 2] = (doublereal)((nbb << 1) + nba);
    }
    /* Test the input arguments */
    if(!notrna && !lsame_(trana, "C", 1, 1))
    {
        *info = -1;
    }
    else if(!notrnb && !lsame_(tranb, "C", 1, 1))
    {
        *info = -2;
    }
    else if(*isgn != 1 && *isgn != -1)
    {
        *info = -3;
    }
    else if(*m < 0)
    {
        *info = -4;
    }
    else if(*n < 0)
    {
        *info = -5;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -7;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -9;
    }
    else if(*ldc < fla_max(1, *m))
    {
        *info = -11;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZTRSYL3", &i__1, (ftnlen)7);
        free(wnrm);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        free(wnrm);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    *scale = 1.;
    if(*m == 0 || *n == 0)
    {
        free(wnrm);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Use unblocked code for small problems or if insufficient */
    /* workspace is provided */
    if(fla_min(nba, nbb) == 1 || *ldswork < fla_max(nba, nbb))
    {
        ztrsyl_(trana, tranb, isgn, m, n, &a[a_offset], lda, &b[b_offset], ldb, &c__[c_offset], ldc,
                scale, info);
        free(wnrm);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Set constants to control overflow */
    smlnum = dlamch_("S");
    bignum = 1. / smlnum;
    /* Set local scaling factors. */
    i__1 = nbb;
    for(l = 1; l <= i__1; ++l)
    {
        i__2 = nba;
        for(k = 1; k <= i__2; ++k)
        {
            swork[k + l * swork_dim1] = 1.;
        }
    }
    /* Fallback scaling factor to prevent flushing of SWORK( K, L ) to zero. */
    /* This scaling is to ensure compatibility with TRSYL and may get flushed. */
    buf = 1.;
    /* Compute upper bounds of blocks of A and B */
    awrk = nbb;
    i__1 = nba;
    for(k = 1; k <= i__1; ++k)
    {
        k1 = (k - 1) * nb + 1;
        /* Computing fla_min */
        i__2 = k * nb;
        k2 = fla_min(i__2, *m) + 1;
        i__2 = nba;
        for(l = k; l <= i__2; ++l)
        {
            l1 = (l - 1) * nb + 1;
            /* Computing fla_min */
            i__3 = l * nb;
            l2 = fla_min(i__3, *m) + 1;
            if(notrna)
            {
                i__3 = k2 - k1;
                i__4 = l2 - l1;
                swork[k + (awrk + l) * swork_dim1]
                    = zlange_("I", &i__3, &i__4, &a[k1 + l1 * a_dim1], lda, wnrm);
            }
            else
            {
                i__3 = k2 - k1;
                i__4 = l2 - l1;
                swork[l + (awrk + k) * swork_dim1]
                    = zlange_("1", &i__3, &i__4, &a[k1 + l1 * a_dim1], lda, wnrm);
            }
        }
    }
    bwrk = nbb + nba;
    i__1 = nbb;
    for(k = 1; k <= i__1; ++k)
    {
        k1 = (k - 1) * nb + 1;
        /* Computing fla_min */
        i__2 = k * nb;
        k2 = fla_min(i__2, *n) + 1;
        i__2 = nbb;
        for(l = k; l <= i__2; ++l)
        {
            l1 = (l - 1) * nb + 1;
            /* Computing fla_min */
            i__3 = l * nb;
            l2 = fla_min(i__3, *n) + 1;
            if(notrnb)
            {
                i__3 = k2 - k1;
                i__4 = l2 - l1;
                swork[k + (bwrk + l) * swork_dim1]
                    = zlange_("I", &i__3, &i__4, &b[k1 + l1 * b_dim1], ldb, wnrm);
            }
            else
            {
                i__3 = k2 - k1;
                i__4 = l2 - l1;
                swork[l + (bwrk + k) * swork_dim1]
                    = zlange_("1", &i__3, &i__4, &b[k1 + l1 * b_dim1], ldb, wnrm);
            }
        }
    }
    sgn = (doublereal)(*isgn);
    z__1.r = sgn;
    z__1.i = 0.; // , expr subst
    csgn.r = z__1.r;
    csgn.i = z__1.i; // , expr subst
    if(notrna && notrnb)
    {
        /* Solve A*X + ISGN*X*B = scale*C. */
        /* The (K,L)th block of X is determined starting from */
        /* bottom-left corner column by column by */
        /* A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */
        /* Where */
        /* M L-1 */
        /* R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]. */
        /* I=K+1 J=1 */
        /* Start loop over block rows (index = K) and block columns (index = L) */
        for(k = nba; k >= 1; --k)
        {
            /* K1: row index of the first row in X( K, L ) */
            /* K2: row index of the first row in X( K+1, L ) */
            /* so the K2 - K1 is the column count of the block X( K, L ) */
            k1 = (k - 1) * nb + 1;
            /* Computing fla_min */
            i__1 = k * nb;
            k2 = fla_min(i__1, *m) + 1;
            i__1 = nbb;
            for(l = 1; l <= i__1; ++l)
            {
                /* L1: column index of the first column in X( K, L ) */
                /* L2: column index of the first column in X( K, L + 1) */
                /* so that L2 - L1 is the row count of the block X( K, L ) */
                l1 = (l - 1) * nb + 1;
                /* Computing fla_min */
                i__2 = l * nb;
                l2 = fla_min(i__2, *n) + 1;
                i__2 = k2 - k1;
                i__3 = l2 - l1;
                ztrsyl_(trana, tranb, isgn, &i__2, &i__3, &a[k1 + k1 * a_dim1], lda,
                        &b[l1 + l1 * b_dim1], ldb, &c__[k1 + l1 * c_dim1], ldc, &scaloc, &iinfo);
                *info = fla_max(*info, iinfo);
                if(scaloc * swork[k + l * swork_dim1] == 0.)
                {
                    if(scaloc == 0.)
                    {
                        /* The magnitude of the largest entry of X(K1:K2-1, L1:L2-1) */
                        /* is larger than the product of BIGNUM**2 and cannot be */
                        /* represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1). */
                        /* Mark the computation as pointless. */
                        buf = 0.;
                    }
                    else
                    {
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        buf *= pow_dd(&c_b18, &d__1);
                    }
                    i__2 = nbb;
                    for(jj = 1; jj <= i__2; ++jj)
                    {
                        i__3 = nba;
                        for(ll = 1; ll <= i__3; ++ll)
                        {
                            /* Bound by BIGNUM to not introduce Inf. The value */
                            /* is irrelevant;
                            corresponding entries of the */
                            /* solution will be flushed in consistency scaling. */
                            /* Computing fla_min */
                            frexp(scaloc, (int *)&temp);
                            d__3 = temp;
                            d__1 = bignum;
                            d__2 = swork[ll + jj * swork_dim1]
                                   / pow_dd(&c_b18, &d__3); // , expr subst
                            swork[ll + jj * swork_dim1] = fla_min(d__1, d__2);
                        }
                    }
                }
                swork[k + l * swork_dim1] = scaloc * swork[k + l * swork_dim1];
                i__2 = k2 - k1;
                i__3 = l2 - l1;
                xnrm = zlange_("I", &i__2, &i__3, &c__[k1 + l1 * c_dim1], ldc, wnrm);
                for(i__ = k - 1; i__ >= 1; --i__)
                {
                    /* C( I, L ) := C( I, L ) - A( I, K ) * C( K, L ) */
                    i1 = (i__ - 1) * nb + 1;
                    /* Computing fla_min */
                    i__2 = i__ * nb;
                    i2 = fla_min(i__2, *m) + 1;
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__2 = i2 - i1;
                    i__3 = l2 - l1;
                    cnrm = zlange_("I", &i__2, &i__3, &c__[i1 + l1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[i__ + l * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1, d__2);
                    cnrm *= scamin / swork[i__ + l * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    anrm = swork[i__ + (awrk + k) * swork_dim1];
                    scaloc = dlarmm_(&anrm, &xnrm, &cnrm);
                    if(scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        buf *= pow_dd(&c_b18, &d__1);
                        i__2 = nbb;
                        for(jj = 1; jj <= i__2; ++jj)
                        {
                            i__3 = nba;
                            for(ll = 1; ll <= i__3; ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *)&temp);
                                d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1]
                                       / pow_dd(&c_b18, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1, d__2);
                            }
                        }
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scamin /= pow_dd(&c_b18, &d__1);
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scaloc /= pow_dd(&c_b18, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to C( I, L ) and C( K, L ). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__2 = l2 - 1;
                        for(jj = l1; jj <= i__2; ++jj)
                        {
                            i__3 = k2 - k1;
                            zdscal_(&i__3, &scal, &c__[k1 + jj * c_dim1], &c__1);
                        }
                    }
                    scal = scamin / swork[i__ + l * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__2 = l2 - 1;
                        for(ll = l1; ll <= i__2; ++ll)
                        {
                            i__3 = i2 - i1;
                            zdscal_(&i__3, &scal, &c__[i1 + ll * c_dim1], &c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[i__ + l * swork_dim1] = scamin * scaloc;
                    i__2 = i2 - i1;
                    i__3 = l2 - l1;
                    i__4 = k2 - k1;
                    z__1.r = -1.;
                    z__1.i = -0.; // , expr subst
                    zgemm_("N", "N", &i__2, &i__3, &i__4, &z__1, &a[i1 + k1 * a_dim1], lda,
                           &c__[k1 + l1 * c_dim1], ldc, &c_b1, &c__[i1 + l1 * c_dim1], ldc);
                }
                i__2 = nbb;
                for(j = l + 1; j <= i__2; ++j)
                {
                    /* C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( L, J ) */
                    j1 = (j - 1) * nb + 1;
                    /* Computing fla_min */
                    i__3 = j * nb;
                    j2 = fla_min(i__3, *n) + 1;
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__3 = k2 - k1;
                    i__4 = j2 - j1;
                    cnrm = zlange_("I", &i__3, &i__4, &c__[k1 + j1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[k + j * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1, d__2);
                    cnrm *= scamin / swork[k + j * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    bnrm = swork[l + (bwrk + j) * swork_dim1];
                    scaloc = dlarmm_(&bnrm, &xnrm, &cnrm);
                    if(scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        buf *= pow_dd(&c_b18, &d__1);
                        i__3 = nbb;
                        for(jj = 1; jj <= i__3; ++jj)
                        {
                            i__4 = nba;
                            for(ll = 1; ll <= i__4; ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *)&temp);
                                d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1]
                                       / pow_dd(&c_b18, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1, d__2);
                            }
                        }
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scamin /= pow_dd(&c_b18, &d__1);
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scaloc /= pow_dd(&c_b18, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to C( K, J ) and C( K, L). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__3 = l2 - 1;
                        for(ll = l1; ll <= i__3; ++ll)
                        {
                            i__4 = k2 - k1;
                            zdscal_(&i__4, &scal, &c__[k1 + ll * c_dim1], &c__1);
                        }
                    }
                    scal = scamin / swork[k + j * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__3 = j2 - 1;
                        for(jj = j1; jj <= i__3; ++jj)
                        {
                            i__4 = k2 - k1;
                            zdscal_(&i__4, &scal, &c__[k1 + jj * c_dim1], &c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[k + j * swork_dim1] = scamin * scaloc;
                    i__3 = k2 - k1;
                    i__4 = j2 - j1;
                    i__5 = l2 - l1;
                    z__1.r = -csgn.r;
                    z__1.i = -csgn.i; // , expr subst
                    zgemm_("N", "N", &i__3, &i__4, &i__5, &z__1, &c__[k1 + l1 * c_dim1], ldc,
                           &b[l1 + j1 * b_dim1], ldb, &c_b1, &c__[k1 + j1 * c_dim1], ldc);
                }
            }
        }
    }
    else if(!notrna && notrnb)
    {
        /* Solve A**H *X + ISGN*X*B = scale*C. */
        /* The (K,L)th block of X is determined starting from */
        /* upper-left corner column by column by */
        /* A(K,K)**H*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */
        /* Where */
        /* K-1 L-1 */
        /* R(K,L) = SUM [A(I,K)**H*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)] */
        /* I=1 J=1 */
        /* Start loop over block rows (index = K) and block columns (index = L) */
        i__1 = nba;
        for(k = 1; k <= i__1; ++k)
        {
            /* K1: row index of the first row in X( K, L ) */
            /* K2: row index of the first row in X( K+1, L ) */
            /* so the K2 - K1 is the column count of the block X( K, L ) */
            k1 = (k - 1) * nb + 1;
            /* Computing fla_min */
            i__2 = k * nb;
            k2 = fla_min(i__2, *m) + 1;
            i__2 = nbb;
            for(l = 1; l <= i__2; ++l)
            {
                /* L1: column index of the first column in X( K, L ) */
                /* L2: column index of the first column in X( K, L + 1) */
                /* so that L2 - L1 is the row count of the block X( K, L ) */
                l1 = (l - 1) * nb + 1;
                /* Computing fla_min */
                i__3 = l * nb;
                l2 = fla_min(i__3, *n) + 1;
                i__3 = k2 - k1;
                i__4 = l2 - l1;
                ztrsyl_(trana, tranb, isgn, &i__3, &i__4, &a[k1 + k1 * a_dim1], lda,
                        &b[l1 + l1 * b_dim1], ldb, &c__[k1 + l1 * c_dim1], ldc, &scaloc, &iinfo);
                *info = fla_max(*info, iinfo);
                if(scaloc * swork[k + l * swork_dim1] == 0.)
                {
                    if(scaloc == 0.)
                    {
                        /* The magnitude of the largest entry of X(K1:K2-1, L1:L2-1) */
                        /* is larger than the product of BIGNUM**2 and cannot be */
                        /* represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1). */
                        /* Mark the computation as pointless. */
                        buf = 0.;
                    }
                    else
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        buf *= pow_dd(&c_b18, &d__1);
                    }
                    i__3 = nbb;
                    for(jj = 1; jj <= i__3; ++jj)
                    {
                        i__4 = nba;
                        for(ll = 1; ll <= i__4; ++ll)
                        {
                            /* Bound by BIGNUM to not introduce Inf. The value */
                            /* is irrelevant;
                            corresponding entries of the */
                            /* solution will be flushed in consistency scaling. */
                            /* Computing fla_min */
                            frexp(scaloc, (int *)&temp);
                            d__3 = temp;
                            d__1 = bignum;
                            d__2 = swork[ll + jj * swork_dim1]
                                   / pow_dd(&c_b18, &d__3); // , expr subst
                            swork[ll + jj * swork_dim1] = fla_min(d__1, d__2);
                        }
                    }
                }
                swork[k + l * swork_dim1] = scaloc * swork[k + l * swork_dim1];
                i__3 = k2 - k1;
                i__4 = l2 - l1;
                xnrm = zlange_("I", &i__3, &i__4, &c__[k1 + l1 * c_dim1], ldc, wnrm);
                i__3 = nba;
                for(i__ = k + 1; i__ <= i__3; ++i__)
                {
                    /* C( I, L ) := C( I, L ) - A( K, I )**H * C( K, L ) */
                    i1 = (i__ - 1) * nb + 1;
                    /* Computing fla_min */
                    i__4 = i__ * nb;
                    i2 = fla_min(i__4, *m) + 1;
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__4 = i2 - i1;
                    i__5 = l2 - l1;
                    cnrm = zlange_("I", &i__4, &i__5, &c__[i1 + l1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[i__ + l * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1, d__2);
                    cnrm *= scamin / swork[i__ + l * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    anrm = swork[i__ + (awrk + k) * swork_dim1];
                    scaloc = dlarmm_(&anrm, &xnrm, &cnrm);
                    if(scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        buf *= pow_dd(&c_b18, &d__1);
                        i__4 = nbb;
                        for(jj = 1; jj <= i__4; ++jj)
                        {
                            i__5 = nba;
                            for(ll = 1; ll <= i__5; ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *)&temp);
                                d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1]
                                       / pow_dd(&c_b18, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1, d__2);
                            }
                        }
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scamin /= pow_dd(&c_b18, &d__1);
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scaloc /= pow_dd(&c_b18, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to to C( I, L ) and C( K, L). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__4 = l2 - 1;
                        for(ll = l1; ll <= i__4; ++ll)
                        {
                            i__5 = k2 - k1;
                            zdscal_(&i__5, &scal, &c__[k1 + ll * c_dim1], &c__1);
                        }
                    }
                    scal = scamin / swork[i__ + l * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__4 = l2 - 1;
                        for(ll = l1; ll <= i__4; ++ll)
                        {
                            i__5 = i2 - i1;
                            zdscal_(&i__5, &scal, &c__[i1 + ll * c_dim1], &c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[i__ + l * swork_dim1] = scamin * scaloc;
                    i__4 = i2 - i1;
                    i__5 = l2 - l1;
                    i__6 = k2 - k1;
                    z__1.r = -1.;
                    z__1.i = -0.; // , expr subst
                    zgemm_("C", "N", &i__4, &i__5, &i__6, &z__1, &a[k1 + i1 * a_dim1], lda,
                           &c__[k1 + l1 * c_dim1], ldc, &c_b1, &c__[i1 + l1 * c_dim1], ldc);
                }
                i__3 = nbb;
                for(j = l + 1; j <= i__3; ++j)
                {
                    /* C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( L, J ) */
                    j1 = (j - 1) * nb + 1;
                    /* Computing fla_min */
                    i__4 = j * nb;
                    j2 = fla_min(i__4, *n) + 1;
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__4 = k2 - k1;
                    i__5 = j2 - j1;
                    cnrm = zlange_("I", &i__4, &i__5, &c__[k1 + j1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[k + j * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1, d__2);
                    cnrm *= scamin / swork[k + j * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    bnrm = swork[l + (bwrk + j) * swork_dim1];
                    scaloc = dlarmm_(&bnrm, &xnrm, &cnrm);
                    if(scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        buf *= pow_dd(&c_b18, &d__1);
                        i__4 = nbb;
                        for(jj = 1; jj <= i__4; ++jj)
                        {
                            i__5 = nba;
                            for(ll = 1; ll <= i__5; ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *)&temp);
                                d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1]
                                       / pow_dd(&c_b18, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1, d__2);
                            }
                        }
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scamin /= pow_dd(&c_b18, &d__1);
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scaloc /= pow_dd(&c_b18, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to to C( K, J ) and C( K, L). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__4 = l2 - 1;
                        for(ll = l1; ll <= i__4; ++ll)
                        {
                            i__5 = k2 - k1;
                            zdscal_(&i__5, &scal, &c__[k1 + ll * c_dim1], &c__1);
                        }
                    }
                    scal = scamin / swork[k + j * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__4 = j2 - 1;
                        for(jj = j1; jj <= i__4; ++jj)
                        {
                            i__5 = k2 - k1;
                            zdscal_(&i__5, &scal, &c__[k1 + jj * c_dim1], &c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[k + j * swork_dim1] = scamin * scaloc;
                    i__4 = k2 - k1;
                    i__5 = j2 - j1;
                    i__6 = l2 - l1;
                    z__1.r = -csgn.r;
                    z__1.i = -csgn.i; // , expr subst
                    zgemm_("N", "N", &i__4, &i__5, &i__6, &z__1, &c__[k1 + l1 * c_dim1], ldc,
                           &b[l1 + j1 * b_dim1], ldb, &c_b1, &c__[k1 + j1 * c_dim1], ldc);
                }
            }
        }
    }
    else if(!notrna && !notrnb)
    {
        /* Solve A**H *X + ISGN*X*B**H = scale*C. */
        /* The (K,L)th block of X is determined starting from */
        /* top-right corner column by column by */
        /* A(K,K)**H*X(K,L) + ISGN*X(K,L)*B(L,L)**H = C(K,L) - R(K,L) */
        /* Where */
        /* K-1 N */
        /* R(K,L) = SUM [A(I,K)**H*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**H]. */
        /* I=1 J=L+1 */
        /* Start loop over block rows (index = K) and block columns (index = L) */
        i__1 = nba;
        for(k = 1; k <= i__1; ++k)
        {
            /* K1: row index of the first row in X( K, L ) */
            /* K2: row index of the first row in X( K+1, L ) */
            /* so the K2 - K1 is the column count of the block X( K, L ) */
            k1 = (k - 1) * nb + 1;
            /* Computing fla_min */
            i__2 = k * nb;
            k2 = fla_min(i__2, *m) + 1;
            for(l = nbb; l >= 1; --l)
            {
                /* L1: column index of the first column in X( K, L ) */
                /* L2: column index of the first column in X( K, L + 1) */
                /* so that L2 - L1 is the row count of the block X( K, L ) */
                l1 = (l - 1) * nb + 1;
                /* Computing fla_min */
                i__2 = l * nb;
                l2 = fla_min(i__2, *n) + 1;
                i__2 = k2 - k1;
                i__3 = l2 - l1;
                ztrsyl_(trana, tranb, isgn, &i__2, &i__3, &a[k1 + k1 * a_dim1], lda,
                        &b[l1 + l1 * b_dim1], ldb, &c__[k1 + l1 * c_dim1], ldc, &scaloc, &iinfo);
                *info = fla_max(*info, iinfo);
                if(scaloc * swork[k + l * swork_dim1] == 0.)
                {
                    if(scaloc == 0.)
                    {
                        /* The magnitude of the largest entry of X(K1:K2-1, L1:L2-1) */
                        /* is larger than the product of BIGNUM**2 and cannot be */
                        /* represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1). */
                        /* Mark the computation as pointless. */
                        buf = 0.;
                    }
                    else
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        buf *= pow_dd(&c_b18, &d__1);
                    }
                    i__2 = nbb;
                    for(jj = 1; jj <= i__2; ++jj)
                    {
                        i__3 = nba;
                        for(ll = 1; ll <= i__3; ++ll)
                        {
                            /* Bound by BIGNUM to not introduce Inf. The value */
                            /* is irrelevant;
                            corresponding entries of the */
                            /* solution will be flushed in consistency scaling. */
                            /* Computing fla_min */
                            frexp(scaloc, (int *)&temp);
                            d__3 = temp;
                            d__1 = bignum;
                            d__2 = swork[ll + jj * swork_dim1]
                                   / pow_dd(&c_b18, &d__3); // , expr subst
                            swork[ll + jj * swork_dim1] = fla_min(d__1, d__2);
                        }
                    }
                }
                swork[k + l * swork_dim1] = scaloc * swork[k + l * swork_dim1];
                i__2 = k2 - k1;
                i__3 = l2 - l1;
                xnrm = zlange_("I", &i__2, &i__3, &c__[k1 + l1 * c_dim1], ldc, wnrm);
                i__2 = nba;
                for(i__ = k + 1; i__ <= i__2; ++i__)
                {
                    /* C( I, L ) := C( I, L ) - A( K, I )**H * C( K, L ) */
                    i1 = (i__ - 1) * nb + 1;
                    /* Computing fla_min */
                    i__3 = i__ * nb;
                    i2 = fla_min(i__3, *m) + 1;
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__3 = i2 - i1;
                    i__4 = l2 - l1;
                    cnrm = zlange_("I", &i__3, &i__4, &c__[i1 + l1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[i__ + l * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1, d__2);
                    cnrm *= scamin / swork[i__ + l * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    anrm = swork[i__ + (awrk + k) * swork_dim1];
                    scaloc = dlarmm_(&anrm, &xnrm, &cnrm);
                    if(scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        buf *= pow_dd(&c_b18, &d__1);
                        i__3 = nbb;
                        for(jj = 1; jj <= i__3; ++jj)
                        {
                            i__4 = nba;
                            for(ll = 1; ll <= i__4; ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *)&temp);
                                d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1]
                                       / pow_dd(&c_b18, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1, d__2);
                            }
                        }
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scamin /= pow_dd(&c_b18, &d__1);
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scaloc /= pow_dd(&c_b18, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to C( I, L ) and C( K, L). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__3 = l2 - 1;
                        for(ll = l1; ll <= i__3; ++ll)
                        {
                            i__4 = k2 - k1;
                            zdscal_(&i__4, &scal, &c__[k1 + ll * c_dim1], &c__1);
                        }
                    }
                    scal = scamin / swork[i__ + l * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__3 = l2 - 1;
                        for(ll = l1; ll <= i__3; ++ll)
                        {
                            i__4 = i2 - i1;
                            zdscal_(&i__4, &scal, &c__[i1 + ll * c_dim1], &c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[i__ + l * swork_dim1] = scamin * scaloc;
                    i__3 = i2 - i1;
                    i__4 = l2 - l1;
                    i__5 = k2 - k1;
                    z__1.r = -1.;
                    z__1.i = -0.; // , expr subst
                    zgemm_("C", "N", &i__3, &i__4, &i__5, &z__1, &a[k1 + i1 * a_dim1], lda,
                           &c__[k1 + l1 * c_dim1], ldc, &c_b1, &c__[i1 + l1 * c_dim1], ldc);
                }
                i__2 = l - 1;
                for(j = 1; j <= i__2; ++j)
                {
                    /* C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( J, L )**H */
                    j1 = (j - 1) * nb + 1;
                    /* Computing fla_min */
                    i__3 = j * nb;
                    j2 = fla_min(i__3, *n) + 1;
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__3 = k2 - k1;
                    i__4 = j2 - j1;
                    cnrm = zlange_("I", &i__3, &i__4, &c__[k1 + j1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[k + j * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1, d__2);
                    cnrm *= scamin / swork[k + j * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    bnrm = swork[l + (bwrk + j) * swork_dim1];
                    scaloc = dlarmm_(&bnrm, &xnrm, &cnrm);
                    if(scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        buf *= pow_dd(&c_b18, &d__1);
                        i__3 = nbb;
                        for(jj = 1; jj <= i__3; ++jj)
                        {
                            i__4 = nba;
                            for(ll = 1; ll <= i__4; ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *)&temp);
                                d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1]
                                       / pow_dd(&c_b18, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1, d__2);
                            }
                        }
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scamin /= pow_dd(&c_b18, &d__1);
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scaloc /= pow_dd(&c_b18, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to C( K, J ) and C( K, L). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__3 = l2 - 1;
                        for(ll = l1; ll <= i__3; ++ll)
                        {
                            i__4 = k2 - k1;
                            zdscal_(&i__4, &scal, &c__[k1 + ll * c_dim1], &c__1);
                        }
                    }
                    scal = scamin / swork[k + j * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__3 = j2 - 1;
                        for(jj = j1; jj <= i__3; ++jj)
                        {
                            i__4 = k2 - k1;
                            zdscal_(&i__4, &scal, &c__[k1 + jj * c_dim1], &c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[k + j * swork_dim1] = scamin * scaloc;
                    i__3 = k2 - k1;
                    i__4 = j2 - j1;
                    i__5 = l2 - l1;
                    z__1.r = -csgn.r;
                    z__1.i = -csgn.i; // , expr subst
                    zgemm_("N", "C", &i__3, &i__4, &i__5, &z__1, &c__[k1 + l1 * c_dim1], ldc,
                           &b[j1 + l1 * b_dim1], ldb, &c_b1, &c__[k1 + j1 * c_dim1], ldc);
                }
            }
        }
    }
    else if(notrna && !notrnb)
    {
        /* Solve A*X + ISGN*X*B**H = scale*C. */
        /* The (K,L)th block of X is determined starting from */
        /* bottom-right corner column by column by */
        /* A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)**H = C(K,L) - R(K,L) */
        /* Where */
        /* M N */
        /* R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**H]. */
        /* I=K+1 J=L+1 */
        /* Start loop over block rows (index = K) and block columns (index = L) */
        for(k = nba; k >= 1; --k)
        {
            /* K1: row index of the first row in X( K, L ) */
            /* K2: row index of the first row in X( K+1, L ) */
            /* so the K2 - K1 is the column count of the block X( K, L ) */
            k1 = (k - 1) * nb + 1;
            /* Computing fla_min */
            i__1 = k * nb;
            k2 = fla_min(i__1, *m) + 1;
            for(l = nbb; l >= 1; --l)
            {
                /* L1: column index of the first column in X( K, L ) */
                /* L2: column index of the first column in X( K, L + 1) */
                /* so that L2 - L1 is the row count of the block X( K, L ) */
                l1 = (l - 1) * nb + 1;
                /* Computing fla_min */
                i__1 = l * nb;
                l2 = fla_min(i__1, *n) + 1;
                i__1 = k2 - k1;
                i__2 = l2 - l1;
                ztrsyl_(trana, tranb, isgn, &i__1, &i__2, &a[k1 + k1 * a_dim1], lda,
                        &b[l1 + l1 * b_dim1], ldb, &c__[k1 + l1 * c_dim1], ldc, &scaloc, &iinfo);
                *info = fla_max(*info, iinfo);
                if(scaloc * swork[k + l * swork_dim1] == 0.)
                {
                    if(scaloc == 0.)
                    {
                        /* The magnitude of the largest entry of X(K1:K2-1, L1:L2-1) */
                        /* is larger than the product of BIGNUM**2 and cannot be */
                        /* represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1). */
                        /* Mark the computation as pointless. */
                        buf = 0.;
                    }
                    else
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        buf *= pow_dd(&c_b18, &d__1);
                    }
                    i__1 = nbb;
                    for(jj = 1; jj <= i__1; ++jj)
                    {
                        i__2 = nba;
                        for(ll = 1; ll <= i__2; ++ll)
                        {
                            /* Bound by BIGNUM to not introduce Inf. The value */
                            /* is irrelevant;
                            corresponding entries of the */
                            /* solution will be flushed in consistency scaling. */
                            /* Computing fla_min */
                            frexp(scaloc, (int *)&temp);
                            d__3 = temp;
                            d__1 = bignum;
                            d__2 = swork[ll + jj * swork_dim1]
                                   / pow_dd(&c_b18, &d__3); // , expr subst
                            swork[ll + jj * swork_dim1] = fla_min(d__1, d__2);
                        }
                    }
                }
                swork[k + l * swork_dim1] = scaloc * swork[k + l * swork_dim1];
                i__1 = k2 - k1;
                i__2 = l2 - l1;
                xnrm = zlange_("I", &i__1, &i__2, &c__[k1 + l1 * c_dim1], ldc, wnrm);
                i__1 = k - 1;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    /* C( I, L ) := C( I, L ) - A( I, K ) * C( K, L ) */
                    i1 = (i__ - 1) * nb + 1;
                    /* Computing fla_min */
                    i__2 = i__ * nb;
                    i2 = fla_min(i__2, *m) + 1;
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__2 = i2 - i1;
                    i__3 = l2 - l1;
                    cnrm = zlange_("I", &i__2, &i__3, &c__[i1 + l1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[i__ + l * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1, d__2);
                    cnrm *= scamin / swork[i__ + l * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    anrm = swork[i__ + (awrk + k) * swork_dim1];
                    scaloc = dlarmm_(&anrm, &xnrm, &cnrm);
                    if(scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        buf *= pow_dd(&c_b18, &d__1);
                        i__2 = nbb;
                        for(jj = 1; jj <= i__2; ++jj)
                        {
                            i__3 = nba;
                            for(ll = 1; ll <= i__3; ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *)&temp);
                                d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1]
                                       / pow_dd(&c_b18, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1, d__2);
                            }
                        }
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scamin /= pow_dd(&c_b18, &d__1);
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scaloc /= pow_dd(&c_b18, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to C( I, L ) and C( K, L). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__2 = l2 - 1;
                        for(ll = l1; ll <= i__2; ++ll)
                        {
                            i__3 = k2 - k1;
                            zdscal_(&i__3, &scal, &c__[k1 + ll * c_dim1], &c__1);
                        }
                    }
                    scal = scamin / swork[i__ + l * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__2 = l2 - 1;
                        for(ll = l1; ll <= i__2; ++ll)
                        {
                            i__3 = i2 - i1;
                            zdscal_(&i__3, &scal, &c__[i1 + ll * c_dim1], &c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[i__ + l * swork_dim1] = scamin * scaloc;
                    i__2 = i2 - i1;
                    i__3 = l2 - l1;
                    i__4 = k2 - k1;
                    z__1.r = -1.;
                    z__1.i = -0.; // , expr subst
                    zgemm_("N", "N", &i__2, &i__3, &i__4, &z__1, &a[i1 + k1 * a_dim1], lda,
                           &c__[k1 + l1 * c_dim1], ldc, &c_b1, &c__[i1 + l1 * c_dim1], ldc);
                }
                i__1 = l - 1;
                for(j = 1; j <= i__1; ++j)
                {
                    /* C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( J, L )**H */
                    j1 = (j - 1) * nb + 1;
                    /* Computing fla_min */
                    i__2 = j * nb;
                    j2 = fla_min(i__2, *n) + 1;
                    /* Compute scaling factor to survive the linear update */
                    /* simulating consistent scaling. */
                    i__2 = k2 - k1;
                    i__3 = j2 - j1;
                    cnrm = zlange_("I", &i__2, &i__3, &c__[k1 + j1 * c_dim1], ldc, wnrm);
                    /* Computing fla_min */
                    d__1 = swork[k + j * swork_dim1];
                    d__2 = swork[k + l * swork_dim1]; // , expr subst
                    scamin = fla_min(d__1, d__2);
                    cnrm *= scamin / swork[k + j * swork_dim1];
                    xnrm *= scamin / swork[k + l * swork_dim1];
                    bnrm = swork[l + (bwrk + j) * swork_dim1];
                    scaloc = dlarmm_(&bnrm, &xnrm, &cnrm);
                    if(scaloc * scamin == 0.)
                    {
                        /* Use second scaling factor to prevent flushing to zero. */
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        buf *= pow_dd(&c_b18, &d__1);
                        i__2 = nbb;
                        for(jj = 1; jj <= i__2; ++jj)
                        {
                            i__3 = nba;
                            for(ll = 1; ll <= i__3; ++ll)
                            {
                                /* Computing fla_min */
                                frexp(scaloc, (int *)&temp);
                                d__3 = temp;
                                d__1 = bignum;
                                d__2 = swork[ll + jj * swork_dim1]
                                       / pow_dd(&c_b18, &d__3); // , expr subst
                                swork[ll + jj * swork_dim1] = fla_min(d__1, d__2);
                            }
                        }
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scamin /= pow_dd(&c_b18, &d__1);
                        frexp(scaloc, (int *)&temp);
                        d__1 = temp;
                        scaloc /= pow_dd(&c_b18, &d__1);
                    }
                    cnrm *= scaloc;
                    xnrm *= scaloc;
                    /* Simultaneously apply the robust update factor and the */
                    /* consistency scaling factor to C( K, J ) and C( K, L). */
                    scal = scamin / swork[k + l * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__2 = l2 - 1;
                        for(jj = l1; jj <= i__2; ++jj)
                        {
                            i__3 = k2 - k1;
                            zdscal_(&i__3, &scal, &c__[k1 + jj * c_dim1], &c__1);
                        }
                    }
                    scal = scamin / swork[k + j * swork_dim1] * scaloc;
                    if(scal != 1.)
                    {
                        i__2 = j2 - 1;
                        for(jj = j1; jj <= i__2; ++jj)
                        {
                            i__3 = k2 - k1;
                            zdscal_(&i__3, &scal, &c__[k1 + jj * c_dim1], &c__1);
                        }
                    }
                    /* Record current scaling factor */
                    swork[k + l * swork_dim1] = scamin * scaloc;
                    swork[k + j * swork_dim1] = scamin * scaloc;
                    i__2 = k2 - k1;
                    i__3 = j2 - j1;
                    i__4 = l2 - l1;
                    z__1.r = -csgn.r;
                    z__1.i = -csgn.i; // , expr subst
                    zgemm_("N", "C", &i__2, &i__3, &i__4, &z__1, &c__[k1 + l1 * c_dim1], ldc,
                           &b[j1 + l1 * b_dim1], ldb, &c_b1, &c__[k1 + j1 * c_dim1], ldc);
                }
            }
        }
    }
    /* Reduce local scaling factors */
    *scale = swork[swork_dim1 + 1];
    i__1 = nba;
    for(k = 1; k <= i__1; ++k)
    {
        i__2 = nbb;
        for(l = 1; l <= i__2; ++l)
        {
            /* Computing fla_min */
            d__1 = *scale;
            d__2 = swork[k + l * swork_dim1]; // , expr subst
            *scale = fla_min(d__1, d__2);
        }
    }
    if(*scale == 0.)
    {
        /* The magnitude of the largest entry of the solution is larger */
        /* than the product of BIGNUM**2 and cannot be represented in the */
        /* form (1/SCALE)*X if SCALE is DOUBLE PRECISION. Set SCALE to */
        /* zero and give up. */
        swork[swork_dim1 + 1] = (doublereal)fla_max(nba, nbb);
        swork[swork_dim1 + 2] = (doublereal)((nbb << 1) + nba);
        free(wnrm);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Realize consistent scaling */
    i__1 = nba;
    for(k = 1; k <= i__1; ++k)
    {
        k1 = (k - 1) * nb + 1;
        /* Computing fla_min */
        i__2 = k * nb;
        k2 = fla_min(i__2, *m) + 1;
        i__2 = nbb;
        for(l = 1; l <= i__2; ++l)
        {
            l1 = (l - 1) * nb + 1;
            /* Computing fla_min */
            i__3 = l * nb;
            l2 = fla_min(i__3, *n) + 1;
            scal = *scale / swork[k + l * swork_dim1];
            if(scal != 1.)
            {
                i__3 = l2 - 1;
                for(ll = l1; ll <= i__3; ++ll)
                {
                    i__4 = k2 - k1;
                    zdscal_(&i__4, &scal, &c__[k1 + ll * c_dim1], &c__1);
                }
            }
        }
    }
    if(buf != 1. && buf > 0.)
    {
        /* Decrease SCALE as much as possible. */
        /* Computing fla_min */
        d__1 = *scale / smlnum;
        d__2 = 1. / buf; // , expr subst
        scaloc = fla_min(d__1, d__2);
        buf *= scaloc;
        *scale /= scaloc;
    }
    if(buf != 1. && buf > 0.)
    {
        /* In case of overly aggressive scaling during the computation, */
        /* flushing of the global scale factor may be prevented by */
        /* undoing some of the scaling. This step is to ensure that */
        /* this routine flushes only scale factors that TRSYL also */
        /* flushes and be usable as a drop-in replacement. */
        /* How much can the normwise largest entry be upscaled? */
        /* Computing fla_max */
        i__1 = c_dim1 + 1;
        d__3 = (d__1 = c__[i__1].r, f2c_abs(d__1));
        d__4 = (d__2 = d_imag(&c__[c_dim1 + 1]), f2c_abs(d__2)); // , expr subst
        scal = fla_max(d__3, d__4);
        i__1 = *m;
        for(k = 1; k <= i__1; ++k)
        {
            i__2 = *n;
            for(l = 1; l <= i__2; ++l)
            {
                /* Computing fla_max */
                i__3 = k + l * c_dim1;
                d__3 = scal, d__4 = (d__1 = c__[i__3].r, f2c_abs(d__1));
                d__3 = fla_max(d__3, d__4);
                d__4 = (d__2 = d_imag(&c__[k + l * c_dim1]), f2c_abs(d__2)); // ; expr subst
                scal = fla_max(d__3, d__4);
            }
        }
        /* Increase BUF as close to 1 as possible and apply scaling. */
        /* Computing fla_min */
        d__1 = bignum / scal;
        d__2 = 1. / buf; // , expr subst
        scaloc = fla_min(d__1, d__2);
        buf *= scaloc;
        zlascl_("G", &c_n1, &c_n1, &c_b106, &scaloc, m, n, &c__[c_offset], ldc, &iinfo);
    }
    /* Combine with buffer scaling factor. SCALE will be flushed if */
    /* BUF is less than one here. */
    *scale *= buf;
    /* Restore workspace dimensions */
    swork[swork_dim1 + 1] = (doublereal)fla_max(nba, nbb);
    swork[swork_dim1 + 2] = (doublereal)((nbb << 1) + nba);
    free(wnrm);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZTRSYL3 */
}
/* ztrsyl3_ */
