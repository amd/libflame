/* ../netlib/v3.9.0/zhetrf_aa.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b2 = {1., 0.};
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b ZHETRF_AA */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHETRF_AA + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrf_
 * aa.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrf_
 * aa.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrf_
 * aa.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHETRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER N, LDA, LWORK, INFO */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHETRF_AA computes the factorization of a complex hermitian matrix A */
/* > using the Aasen's algorithm. The form of the factorization is */
/* > */
/* > A = U**H*T*U or A = L*T*L**H */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and T is a hermitian tridiagonal matrix. */
/* > */
/* > This is the blocked version of the algorithm, calling Level 3 BLAS. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the hermitian matrix A. If UPLO = 'U', the leading */
/* > N-by-N upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading N-by-N lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > */
/* > On exit, the tridiagonal matrix is stored in the diagonals */
/* > and the subdiagonals of A just below (or above) the diagonals, */
/* > and L is stored below (or above) the subdiaonals, when UPLO */
/* > is 'L' (or 'U'). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > On exit, it contains the details of the interchanges, i.e., */
/* > the row and column k of A were interchanged with the */
/* > row and column IPIV(k). */
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
/* > The length of WORK. LWORK >= MAX(1,2*N). For optimum performance */
/* > LWORK >= N*(1+NB), where NB is the optimal blocksize. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2017 */
/* > \ingroup complex16HEcomputational */
/* ===================================================================== */
/* Subroutine */
void zhetrf_aa_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *ipiv,
                doublecomplex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhetrf_aa inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS ", lwork %" FLA_IS
                      "",
                      *uplo, *n, *lda, *lwork);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    integer j;
    extern /* Subroutine */
        void
        zlahef_aa_(char *, integer *, integer *, integer *, doublecomplex *, integer *, integer *,
                   doublecomplex *, integer *, doublecomplex *);
    integer j1, k1, k2, j2, j3, jb, nb, mj, nj;
    doublecomplex alpha;
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        zscal_(integer *, doublecomplex *, doublecomplex *, integer *),
        zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *, doublecomplex *,
               integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *);
    logical upper;
    extern /* Subroutine */
        void
        zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *),
        zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer lwkopt;
    logical lquery;
    /* -- LAPACK computational routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2017 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Determine the block size */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --work;
    /* Function Body */
    nb = ilaenv_(&c__1, "ZHETRF_AA", uplo, n, &c_n1, &c_n1, &c_n1);
    /* Test the input parameters. */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    lquery = *lwork == -1;
    if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -4;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n << 1; // , expr subst
        if(*lwork < fla_max(i__1, i__2) && !lquery)
        {
            *info = -7;
        }
    }
    if(*info == 0)
    {
        lwkopt = (nb + 1) * *n;
        work[1].r = (doublereal)lwkopt;
        work[1].i = 0.; // , expr subst
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHETRF_AA", &i__1, (ftnlen)9);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    ipiv[1] = 1;
    if(*n == 1)
    {
        i__1 = a_dim1 + 1;
        i__2 = a_dim1 + 1;
        d__1 = a[i__2].r;
        a[i__1].r = d__1;
        a[i__1].i = 0.; // , expr subst
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Adjust block size based on the workspace size */
    if(*lwork < (nb + 1) * *n)
    {
        nb = (*lwork - *n) / *n;
    }
    if(upper)
    {
        /* ..................................................... */
        /* Factorize A as U**H*D*U using the upper triangle of A */
        /* ..................................................... */
        /* copy first row A(1, 1:N) into H(1:n) (stored in WORK(1:N)) */
        zcopy_(n, &a[a_dim1 + 1], lda, &work[1], &c__1);
        /* J is the main loop index, increasing from 1 to N in steps of */
        /* JB, where JB is the number of columns factorized by ZLAHEF;
         */
        /* JB is either NB, or N-J+1 for the last block */
        j = 0;
    L10:
        if(j >= *n)
        {
            goto L20;
        }
        /* each step of the main loop */
        /* J is the last column of the previous panel */
        /* J1 is the first column of the current panel */
        /* K1 identifies if the previous column of the panel has been */
        /* explicitly stored, e.g., K1=1 for the first panel, and */
        /* K1=0 for the rest */
        j1 = j + 1;
        /* Computing MIN */
        i__1 = *n - j1 + 1;
        jb = fla_min(i__1, nb);
        k1 = fla_max(1, j) - j;
        /* Panel factorization */
        i__1 = 2 - k1;
        i__2 = *n - j;
        zlahef_aa_(uplo, &i__1, &i__2, &jb, &a[fla_max(1, j) + (j + 1) * a_dim1], lda, &ipiv[j + 1],
                   &work[1], n, &work[*n * nb + 1]);
        /* Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot) */
        /* Computing MIN */
        i__2 = *n;
        i__3 = j + jb + 1; // , expr subst
        i__1 = fla_min(i__2, i__3);
        for(j2 = j + 2; j2 <= i__1; ++j2)
        {
            ipiv[j2] += j;
            if(j2 != ipiv[j2] && j1 - k1 > 2)
            {
                i__2 = j1 - k1 - 2;
                zswap_(&i__2, &a[j2 * a_dim1 + 1], &c__1, &a[ipiv[j2] * a_dim1 + 1], &c__1);
            }
        }
        j += jb;
        /* Trailing submatrix update, where */
        /* the row A(J1-1, J2-1:N) stores U(J1, J2+1:N) and */
        /* WORK stores the current block of the auxiriarly matrix H */
        if(j < *n)
        {
            /* if the first panel and JB=1 (NB=1), then nothing to do */
            if(j1 > 1 || jb > 1)
            {
                /* Merge rank-1 update with BLAS-3 update */
                d_cnjg(&z__1, &a[j + (j + 1) * a_dim1]);
                alpha.r = z__1.r;
                alpha.i = z__1.i; // , expr subst
                i__1 = j + (j + 1) * a_dim1;
                a[i__1].r = 1.;
                a[i__1].i = 0.; // , expr subst
                i__1 = *n - j;
                zcopy_(&i__1, &a[j - 1 + (j + 1) * a_dim1], lda, &work[j + 1 - j1 + 1 + jb * *n],
                       &c__1);
                i__1 = *n - j;
                zscal_(&i__1, &alpha, &work[j + 1 - j1 + 1 + jb * *n], &c__1);
                /* K1 identifies if the previous column of the panel has been */
                /* explicitly stored, e.g., K1=0 and K2=1 for the first panel, */
                /* and K1=1 and K2=0 for the rest */
                if(j1 > 1)
                {
                    /* Not first panel */
                    k2 = 1;
                }
                else
                {
                    /* First panel */
                    k2 = 0;
                    /* First update skips the first column */
                    --jb;
                }
                i__1 = *n;
                i__2 = nb;
                for(j2 = j + 1; i__2 < 0 ? j2 >= i__1 : j2 <= i__1; j2 += i__2)
                {
                    /* Computing MIN */
                    i__3 = nb;
                    i__4 = *n - j2 + 1; // , expr subst
                    nj = fla_min(i__3, i__4);
                    /* Update (J2, J2) diagonal block with ZGEMV */
                    j3 = j2;
                    for(mj = nj - 1; mj >= 1; --mj)
                    {
                        i__3 = jb + 1;
                        z__1.r = -1.;
                        z__1.i = -0.; // , expr subst
                        zgemm_("Conjugate transpose", "Transpose", &c__1, &mj, &i__3, &z__1,
                               &a[j1 - k2 + j3 * a_dim1], lda, &work[j3 - j1 + 1 + k1 * *n], n,
                               &c_b2, &a[j3 + j3 * a_dim1], lda);
                        ++j3;
                    }
                    /* Update off-diagonal block of J2-th block row with ZGEMM */
                    i__3 = *n - j3 + 1;
                    i__4 = jb + 1;
                    z__1.r = -1.;
                    z__1.i = -0.; // , expr subst
                    zgemm_("Conjugate transpose", "Transpose", &nj, &i__3, &i__4, &z__1,
                           &a[j1 - k2 + j2 * a_dim1], lda, &work[j3 - j1 + 1 + k1 * *n], n, &c_b2,
                           &a[j2 + j3 * a_dim1], lda);
                }
                /* Recover T( J, J+1 ) */
                i__2 = j + (j + 1) * a_dim1;
                d_cnjg(&z__1, &alpha);
                a[i__2].r = z__1.r;
                a[i__2].i = z__1.i; // , expr subst
            }
            /* WORK(J+1, 1) stores H(J+1, 1) */
            i__2 = *n - j;
            zcopy_(&i__2, &a[j + 1 + (j + 1) * a_dim1], lda, &work[1], &c__1);
        }
        goto L10;
    }
    else
    {
        /* ..................................................... */
        /* Factorize A as L*D*L**H using the lower triangle of A */
        /* ..................................................... */
        /* copy first column A(1:N, 1) into H(1:N, 1) */
        /* (stored in WORK(1:N)) */
        zcopy_(n, &a[a_dim1 + 1], &c__1, &work[1], &c__1);
        /* J is the main loop index, increasing from 1 to N in steps of */
        /* JB, where JB is the number of columns factorized by ZLAHEF;
         */
        /* JB is either NB, or N-J+1 for the last block */
        j = 0;
    L11:
        if(j >= *n)
        {
            goto L20;
        }
        /* each step of the main loop */
        /* J is the last column of the previous panel */
        /* J1 is the first column of the current panel */
        /* K1 identifies if the previous column of the panel has been */
        /* explicitly stored, e.g., K1=1 for the first panel, and */
        /* K1=0 for the rest */
        j1 = j + 1;
        /* Computing MIN */
        i__2 = *n - j1 + 1;
        jb = fla_min(i__2, nb);
        k1 = fla_max(1, j) - j;
        /* Panel factorization */
        i__2 = 2 - k1;
        i__1 = *n - j;
        zlahef_aa_(uplo, &i__2, &i__1, &jb, &a[j + 1 + fla_max(1, j) * a_dim1], lda, &ipiv[j + 1],
                   &work[1], n, &work[*n * nb + 1]);
        /* Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot) */
        /* Computing MIN */
        i__1 = *n;
        i__3 = j + jb + 1; // , expr subst
        i__2 = fla_min(i__1, i__3);
        for(j2 = j + 2; j2 <= i__2; ++j2)
        {
            ipiv[j2] += j;
            if(j2 != ipiv[j2] && j1 - k1 > 2)
            {
                i__1 = j1 - k1 - 2;
                zswap_(&i__1, &a[j2 + a_dim1], lda, &a[ipiv[j2] + a_dim1], lda);
            }
        }
        j += jb;
        /* Trailing submatrix update, where */
        /* A(J2+1, J1-1) stores L(J2+1, J1) and */
        /* WORK(J2+1, 1) stores H(J2+1, 1) */
        if(j < *n)
        {
            /* if the first panel and JB=1 (NB=1), then nothing to do */
            if(j1 > 1 || jb > 1)
            {
                /* Merge rank-1 update with BLAS-3 update */
                d_cnjg(&z__1, &a[j + 1 + j * a_dim1]);
                alpha.r = z__1.r;
                alpha.i = z__1.i; // , expr subst
                i__2 = j + 1 + j * a_dim1;
                a[i__2].r = 1.;
                a[i__2].i = 0.; // , expr subst
                i__2 = *n - j;
                zcopy_(&i__2, &a[j + 1 + (j - 1) * a_dim1], &c__1, &work[j + 1 - j1 + 1 + jb * *n],
                       &c__1);
                i__2 = *n - j;
                zscal_(&i__2, &alpha, &work[j + 1 - j1 + 1 + jb * *n], &c__1);
                /* K1 identifies if the previous column of the panel has been */
                /* explicitly stored, e.g., K1=0 and K2=1 for the first panel, */
                /* and K1=1 and K2=0 for the rest */
                if(j1 > 1)
                {
                    /* Not first panel */
                    k2 = 1;
                }
                else
                {
                    /* First panel */
                    k2 = 0;
                    /* First update skips the first column */
                    --jb;
                }
                i__2 = *n;
                i__1 = nb;
                for(j2 = j + 1; i__1 < 0 ? j2 >= i__2 : j2 <= i__2; j2 += i__1)
                {
                    /* Computing MIN */
                    i__3 = nb;
                    i__4 = *n - j2 + 1; // , expr subst
                    nj = fla_min(i__3, i__4);
                    /* Update (J2, J2) diagonal block with ZGEMV */
                    j3 = j2;
                    for(mj = nj - 1; mj >= 1; --mj)
                    {
                        i__3 = jb + 1;
                        z__1.r = -1.;
                        z__1.i = -0.; // , expr subst
                        zgemm_("No transpose", "Conjugate transpose", &mj, &c__1, &i__3, &z__1,
                               &work[j3 - j1 + 1 + k1 * *n], n, &a[j3 + (j1 - k2) * a_dim1], lda,
                               &c_b2, &a[j3 + j3 * a_dim1], lda);
                        ++j3;
                    }
                    /* Update off-diagonal block of J2-th block column with ZGEMM */
                    i__3 = *n - j3 + 1;
                    i__4 = jb + 1;
                    z__1.r = -1.;
                    z__1.i = -0.; // , expr subst
                    zgemm_("No transpose", "Conjugate transpose", &i__3, &nj, &i__4, &z__1,
                           &work[j3 - j1 + 1 + k1 * *n], n, &a[j2 + (j1 - k2) * a_dim1], lda, &c_b2,
                           &a[j3 + j2 * a_dim1], lda);
                }
                /* Recover T( J+1, J ) */
                i__1 = j + 1 + j * a_dim1;
                d_cnjg(&z__1, &alpha);
                a[i__1].r = z__1.r;
                a[i__1].i = z__1.i; // , expr subst
            }
            /* WORK(J+1, 1) stores H(J+1, 1) */
            i__1 = *n - j;
            zcopy_(&i__1, &a[j + 1 + (j + 1) * a_dim1], &c__1, &work[1], &c__1);
        }
        goto L11;
    }
L20:
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZHETRF_AA */
}
/* zhetrf_aa__ */
