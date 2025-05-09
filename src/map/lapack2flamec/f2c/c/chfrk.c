/* ../netlib/chfrk.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CHFRK performs a Hermitian rank-k operation for matrix in RFP format. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHFRK + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chfrk.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chfrk.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chfrk.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, */
/* C ) */
/* .. Scalar Arguments .. */
/* REAL ALPHA, BETA */
/* INTEGER K, LDA, N */
/* CHARACTER TRANS, TRANSR, UPLO */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), C( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > Level 3 BLAS like routine for C in RFP Format. */
/* > */
/* > CHFRK performs one of the Hermitian rank--k operations */
/* > */
/* > C := alpha*A*A**H + beta*C, */
/* > */
/* > or */
/* > */
/* > C := alpha*A**H*A + beta*C, */
/* > */
/* > where alpha and beta are real scalars, C is an n--by--n Hermitian */
/* > matrix and A is an n--by--k matrix in the first case and a k--by--n */
/* > matrix in the second case. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANSR */
/* > \verbatim */
/* > TRANSR is CHARACTER*1 */
/* > = 'N': The Normal Form of RFP A is stored;
 */
/* > = 'C': The Conjugate-transpose Form of RFP A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > On entry, UPLO specifies whether the upper or lower */
/* > triangular part of the array C is to be referenced as */
/* > follows: */
/* > */
/* > UPLO = 'U' or 'u' Only the upper triangular part of C */
/* > is to be referenced. */
/* > */
/* > UPLO = 'L' or 'l' Only the lower triangular part of C */
/* > is to be referenced. */
/* > */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > On entry, TRANS specifies the operation to be performed as */
/* > follows: */
/* > */
/* > TRANS = 'N' or 'n' C := alpha*A*A**H + beta*C. */
/* > */
/* > TRANS = 'C' or 'c' C := alpha*A**H*A + beta*C. */
/* > */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > On entry, N specifies the order of the matrix C. N must be */
/* > at least zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > On entry with TRANS = 'N' or 'n', K specifies the number */
/* > of columns of the matrix A, and on entry with */
/* > TRANS = 'C' or 'c', K specifies the number of rows of the */
/* > matrix A. K must be at least zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* > ALPHA is REAL */
/* > On entry, ALPHA specifies the scalar alpha. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,ka) */
/* > where KA */
/* > is K when TRANS = 'N' or 'n', and is N otherwise. Before */
/* > entry with TRANS = 'N' or 'n', the leading N--by--K part of */
/* > the array A must contain the matrix A, otherwise the leading */
/* > K--by--N part of the array A must contain the matrix A. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > On entry, LDA specifies the first dimension of A as declared */
/* > in the calling (sub) program. When TRANS = 'N' or 'n' */
/* > then LDA must be at least fla_max( 1, n ), otherwise LDA must */
/* > be at least fla_max( 1, k ). */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* > BETA is REAL */
/* > On entry, BETA specifies the scalar beta. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX array, dimension (N*(N+1)/2) */
/* > On entry, the matrix A in RFP Format. RFP Format is */
/* > described by TRANSR, UPLO and N. Note that the imaginary */
/* > parts of the diagonal elements need not be set, they are */
/* > assumed to be zero, and on exit they are set to zero. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
void chfrk_(char *transr, char *uplo, char *trans, integer *n, integer *k, real *alpha, complex *a,
            integer *lda, real *beta, complex *c__)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "chfrk inputs: transr %c, uplo %c, trans %c, n %lld, k %lld, lda %lld",
             *transr, *uplo, *trans, *n, *k, *lda);
#else
    snprintf(buffer, 256, "chfrk inputs: transr %c, uplo %c, trans %c, n %d, k %d, lda %d", *transr,
             *uplo, *trans, *n, *k, *lda);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    complex q__1;
    /* Local variables */
    integer j, n1, n2, nk, info;
    complex cbeta;
    logical normaltransr;
    extern /* Subroutine */
        void
        cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *,
               complex *, integer *, complex *, complex *, integer *),
        cherk_(char *, char *, integer *, integer *, real *, complex *, integer *, real *,
               complex *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    integer nrowa;
    logical lower;
    complex calpha;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical nisodd, notrans;
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. */
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --c__;
    /* Function Body */
    info = 0;
    normaltransr = lsame_(transr, "N", 1, 1);
    lower = lsame_(uplo, "L", 1, 1);
    notrans = lsame_(trans, "N", 1, 1);
    if(notrans)
    {
        nrowa = *n;
    }
    else
    {
        nrowa = *k;
    }
    if(!normaltransr && !lsame_(transr, "C", 1, 1))
    {
        info = -1;
    }
    else if(!lower && !lsame_(uplo, "U", 1, 1))
    {
        info = -2;
    }
    else if(!notrans && !lsame_(trans, "C", 1, 1))
    {
        info = -3;
    }
    else if(*n < 0)
    {
        info = -4;
    }
    else if(*k < 0)
    {
        info = -5;
    }
    else if(*lda < fla_max(1, nrowa))
    {
        info = -8;
    }
    if(info != 0)
    {
        i__1 = -info;
        xerbla_("CHFRK ", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible. */
    /* The quick return case: ((ALPHA.EQ.0).AND.(BETA.NE.ZERO)) is not */
    /* done (it is in CHERK for example) and left in the general case. */
    if(*n == 0 || (*alpha == 0.f || *k == 0) && *beta == 1.f)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(*alpha == 0.f && *beta == 0.f)
    {
        i__1 = *n * (*n + 1) / 2;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j;
            c__[i__2].r = 0.f;
            c__[i__2].i = 0.f; // , expr subst
        }
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    q__1.r = *alpha;
    q__1.i = 0.f; // , expr subst
    calpha.r = q__1.r;
    calpha.i = q__1.i; // , expr subst
    q__1.r = *beta;
    q__1.i = 0.f; // , expr subst
    cbeta.r = q__1.r;
    cbeta.i = q__1.i; // , expr subst
    /* C is N-by-N. */
    /* If N is odd, set NISODD = .TRUE., and N1 and N2. */
    /* If N is even, NISODD = .FALSE., and NK. */
    if(*n % 2 == 0)
    {
        nisodd = FALSE_;
        nk = *n / 2;
    }
    else
    {
        nisodd = TRUE_;
        if(lower)
        {
            n2 = *n / 2;
            n1 = *n - n2;
        }
        else
        {
            n1 = *n / 2;
            n2 = *n - n1;
        }
    }
    if(nisodd)
    {
        /* N is odd */
        if(normaltransr)
        {
            /* N is odd and TRANSR = 'N' */
            if(lower)
            {
                /* N is odd, TRANSR = 'N', and UPLO = 'L' */
                if(notrans)
                {
                    /* N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'N' */
                    cherk_("L", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[1], n);
                    cherk_("U", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, beta, &c__[*n + 1],
                           n);
                    cgemm_("N", "C", &n2, &n1, k, &calpha, &a[n1 + 1 + a_dim1], lda, &a[a_dim1 + 1],
                           lda, &cbeta, &c__[n1 + 1], n);
                }
                else
                {
                    /* N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'C' */
                    cherk_("L", "C", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[1], n);
                    cherk_("U", "C", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1], lda, beta,
                           &c__[*n + 1], n);
                    cgemm_("C", "N", &n2, &n1, k, &calpha, &a[(n1 + 1) * a_dim1 + 1], lda,
                           &a[a_dim1 + 1], lda, &cbeta, &c__[n1 + 1], n);
                }
            }
            else
            {
                /* N is odd, TRANSR = 'N', and UPLO = 'U' */
                if(notrans)
                {
                    /* N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'N' */
                    cherk_("L", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[n2 + 1], n);
                    cherk_("U", "N", &n2, k, alpha, &a[n2 + a_dim1], lda, beta, &c__[n1 + 1], n);
                    cgemm_("N", "C", &n1, &n2, k, &calpha, &a[a_dim1 + 1], lda, &a[n2 + a_dim1],
                           lda, &cbeta, &c__[1], n);
                }
                else
                {
                    /* N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'C' */
                    cherk_("L", "C", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[n2 + 1], n);
                    cherk_("U", "C", &n2, k, alpha, &a[n2 * a_dim1 + 1], lda, beta, &c__[n1 + 1],
                           n);
                    cgemm_("C", "N", &n1, &n2, k, &calpha, &a[a_dim1 + 1], lda, &a[n2 * a_dim1 + 1],
                           lda, &cbeta, &c__[1], n);
                }
            }
        }
        else
        {
            /* N is odd, and TRANSR = 'C' */
            if(lower)
            {
                /* N is odd, TRANSR = 'C', and UPLO = 'L' */
                if(notrans)
                {
                    /* N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'N' */
                    cherk_("U", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[1], &n1);
                    cherk_("L", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, beta, &c__[2], &n1);
                    cgemm_("N", "C", &n1, &n2, k, &calpha, &a[a_dim1 + 1], lda, &a[n1 + 1 + a_dim1],
                           lda, &cbeta, &c__[n1 * n1 + 1], &n1);
                }
                else
                {
                    /* N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'C' */
                    cherk_("U", "C", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[1], &n1);
                    cherk_("L", "C", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1], lda, beta, &c__[2],
                           &n1);
                    cgemm_("C", "N", &n1, &n2, k, &calpha, &a[a_dim1 + 1], lda,
                           &a[(n1 + 1) * a_dim1 + 1], lda, &cbeta, &c__[n1 * n1 + 1], &n1);
                }
            }
            else
            {
                /* N is odd, TRANSR = 'C', and UPLO = 'U' */
                if(notrans)
                {
                    /* N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'N' */
                    cherk_("U", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[n2 * n2 + 1],
                           &n2);
                    cherk_("L", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, beta,
                           &c__[n1 * n2 + 1], &n2);
                    cgemm_("N", "C", &n2, &n1, k, &calpha, &a[n1 + 1 + a_dim1], lda, &a[a_dim1 + 1],
                           lda, &cbeta, &c__[1], &n2);
                }
                else
                {
                    /* N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'C' */
                    cherk_("U", "C", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[n2 * n2 + 1],
                           &n2);
                    cherk_("L", "C", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1], lda, beta,
                           &c__[n1 * n2 + 1], &n2);
                    cgemm_("C", "N", &n2, &n1, k, &calpha, &a[(n1 + 1) * a_dim1 + 1], lda,
                           &a[a_dim1 + 1], lda, &cbeta, &c__[1], &n2);
                }
            }
        }
    }
    else
    {
        /* N is even */
        if(normaltransr)
        {
            /* N is even and TRANSR = 'N' */
            if(lower)
            {
                /* N is even, TRANSR = 'N', and UPLO = 'L' */
                if(notrans)
                {
                    /* N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'N' */
                    i__1 = *n + 1;
                    cherk_("L", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[2], &i__1);
                    i__1 = *n + 1;
                    cherk_("U", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, beta, &c__[1], &i__1);
                    i__1 = *n + 1;
                    cgemm_("N", "C", &nk, &nk, k, &calpha, &a[nk + 1 + a_dim1], lda, &a[a_dim1 + 1],
                           lda, &cbeta, &c__[nk + 2], &i__1);
                }
                else
                {
                    /* N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'C' */
                    i__1 = *n + 1;
                    cherk_("L", "C", &nk, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[2], &i__1);
                    i__1 = *n + 1;
                    cherk_("U", "C", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1], lda, beta, &c__[1],
                           &i__1);
                    i__1 = *n + 1;
                    cgemm_("C", "N", &nk, &nk, k, &calpha, &a[(nk + 1) * a_dim1 + 1], lda,
                           &a[a_dim1 + 1], lda, &cbeta, &c__[nk + 2], &i__1);
                }
            }
            else
            {
                /* N is even, TRANSR = 'N', and UPLO = 'U' */
                if(notrans)
                {
                    /* N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'N' */
                    i__1 = *n + 1;
                    cherk_("L", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[nk + 2], &i__1);
                    i__1 = *n + 1;
                    cherk_("U", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, beta, &c__[nk + 1],
                           &i__1);
                    i__1 = *n + 1;
                    cgemm_("N", "C", &nk, &nk, k, &calpha, &a[a_dim1 + 1], lda, &a[nk + 1 + a_dim1],
                           lda, &cbeta, &c__[1], &i__1);
                }
                else
                {
                    /* N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'C' */
                    i__1 = *n + 1;
                    cherk_("L", "C", &nk, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[nk + 2], &i__1);
                    i__1 = *n + 1;
                    cherk_("U", "C", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1], lda, beta,
                           &c__[nk + 1], &i__1);
                    i__1 = *n + 1;
                    cgemm_("C", "N", &nk, &nk, k, &calpha, &a[a_dim1 + 1], lda,
                           &a[(nk + 1) * a_dim1 + 1], lda, &cbeta, &c__[1], &i__1);
                }
            }
        }
        else
        {
            /* N is even, and TRANSR = 'C' */
            if(lower)
            {
                /* N is even, TRANSR = 'C', and UPLO = 'L' */
                if(notrans)
                {
                    /* N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'N' */
                    cherk_("U", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[nk + 1], &nk);
                    cherk_("L", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, beta, &c__[1], &nk);
                    cgemm_("N", "C", &nk, &nk, k, &calpha, &a[a_dim1 + 1], lda, &a[nk + 1 + a_dim1],
                           lda, &cbeta, &c__[(nk + 1) * nk + 1], &nk);
                }
                else
                {
                    /* N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'C' */
                    cherk_("U", "C", &nk, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[nk + 1], &nk);
                    cherk_("L", "C", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1], lda, beta, &c__[1],
                           &nk);
                    cgemm_("C", "N", &nk, &nk, k, &calpha, &a[a_dim1 + 1], lda,
                           &a[(nk + 1) * a_dim1 + 1], lda, &cbeta, &c__[(nk + 1) * nk + 1], &nk);
                }
            }
            else
            {
                /* N is even, TRANSR = 'C', and UPLO = 'U' */
                if(notrans)
                {
                    /* N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'N' */
                    cherk_("U", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
                           &c__[nk * (nk + 1) + 1], &nk);
                    cherk_("L", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, beta,
                           &c__[nk * nk + 1], &nk);
                    cgemm_("N", "C", &nk, &nk, k, &calpha, &a[nk + 1 + a_dim1], lda, &a[a_dim1 + 1],
                           lda, &cbeta, &c__[1], &nk);
                }
                else
                {
                    /* N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'C' */
                    cherk_("U", "C", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
                           &c__[nk * (nk + 1) + 1], &nk);
                    cherk_("L", "C", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1], lda, beta,
                           &c__[nk * nk + 1], &nk);
                    cgemm_("C", "N", &nk, &nk, k, &calpha, &a[(nk + 1) * a_dim1 + 1], lda,
                           &a[a_dim1 + 1], lda, &cbeta, &c__[1], &nk);
                }
            }
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CHFRK */
}
/* chfrk_ */
