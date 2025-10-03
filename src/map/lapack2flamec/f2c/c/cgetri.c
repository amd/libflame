/* ./cgetri.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b2 = {1.f, 0.f};
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
static aocl_int64_t c__2 = 2;
/* > \brief \b CGETRI */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGETRI + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgetri.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgetri.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgetri.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGETRI computes the inverse of a matrix using the LU factorization */
/* > computed by CGETRF. */
/* > */
/* > This method inverts U and then computes inv(A) by solving the system */
/* > inv(A)*L = inv(U) for inv(A). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the factors L and U from the factorization */
/* > A = P*L*U as computed by CGETRF. */
/* > On exit, if INFO = 0, the inverse of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > The pivot indices from CGETRF;
for 1<=i<=N, row i of the */
/* > matrix was interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO=0, then WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= fla_max(1,N). */
/* > For optimal performance LWORK >= N*NB, where NB is */
/* > the optimal blocksize returned by ILAENV. */
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
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, U(i,i) is exactly zero;
the matrix is */
/* > singular and its inverse could not be computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup getri */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void cgetri_(aocl_int_t *n, scomplex *a, aocl_int_t *lda, aocl_int_t *ipiv, scomplex *work,
             aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgetri(n, a, lda, ipiv, work, lwork, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cgetri(&n_64, a, &lda_64, ipiv, work, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_cgetri(aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                        scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "cgetri inputs: n %lld, lda %lld, lwork %lld", *n, *lda, *lwork);
#else
    snprintf(buffer, 256, "cgetri inputs: n %d, lda %d, lwork %d", *n, *lda, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1;
    scomplex q__1;
    /* Local variables */
    aocl_int64_t i__, j, jb, nb, jj, jp, nn, iws;
    aocl_int64_t nbmin;
    aocl_int64_t ldwork;
    aocl_int64_t lwkopt;
    logical lquery;
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --work;
    /* Function Body */
    *info = 0;
    nb = aocl_lapack_ilaenv(&c__1, "CGETRI", " ", n, &c_n1, &c_n1, &c_n1);
    lwkopt = *n * nb;
    r__1 = aocl_lapack_sroundup_lwork(&lwkopt);
    work[1].real = r__1;
    work[1].imag = 0.f; // , expr subst
    lquery = *lwork == -1;
    if(*n < 0)
    {
        *info = -1;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -3;
    }
    else if(*lwork < fla_max(1, *n) && !lquery)
    {
        *info = -6;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CGETRI", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Form inv(U). If INFO > 0 from CTRTRI, then U is singular, */
    /* and the inverse is not computed. */
    aocl_lapack_ctrtri("Upper", "Non-unit", n, &a[a_offset], lda, info);
    if(*info > 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    nbmin = 2;
    ldwork = *n;
    if(nb > 1 && nb < *n)
    {
        /* Computing MAX */
        i__1 = ldwork * nb;
        iws = fla_max(i__1, 1);
        if(*lwork < iws)
        {
            nb = *lwork / ldwork;
            /* Computing MAX */
            i__1 = 2;
            i__2 = aocl_lapack_ilaenv(&c__2, "CGETRI", " ", n, &c_n1, &c_n1, &c_n1); // , expr subst
            nbmin = fla_max(i__1, i__2);
        }
    }
    else
    {
        iws = *n;
    }
    /* Solve the equation inv(A)*L = inv(U) for inv(A). */
    if(nb < nbmin || nb >= *n)
    {
        /* Use unblocked code. */
        for(j = *n; j >= 1; --j)
        {
            /* Copy current column of L to WORK and replace with zeros. */
            i__1 = *n;
            for(i__ = j + 1; i__ <= i__1; ++i__)
            {
                i__2 = i__;
                i__3 = i__ + j * a_dim1;
                work[i__2].real = a[i__3].real;
                work[i__2].imag = a[i__3].imag; // , expr subst
                i__2 = i__ + j * a_dim1;
                a[i__2].real = 0.f;
                a[i__2].imag = 0.f; // , expr subst
                /* L10: */
            }
            /* Compute current column of inv(A). */
            if(j < *n)
            {
                i__1 = *n - j;
                q__1.real = -1.f;
                q__1.imag = -0.f; // , expr subst
                aocl_blas_cgemv("No transpose", n, &i__1, &q__1, &a[(j + 1) * a_dim1 + 1], lda,
                                &work[j + 1], &c__1, &c_b2, &a[j * a_dim1 + 1], &c__1);
            }
            /* L20: */
        }
    }
    else
    {
        /* Use blocked code. */
        nn = (*n - 1) / nb * nb + 1;
        i__1 = -nb;
        for(j = nn; i__1 < 0 ? j >= 1 : j <= 1; j += i__1)
        {
            /* Computing MIN */
            i__2 = nb;
            i__3 = *n - j + 1; // , expr subst
            jb = fla_min(i__2, i__3);
            /* Copy current block column of L to WORK and replace with */
            /* zeros. */
            i__2 = j + jb - 1;
            for(jj = j; jj <= i__2; ++jj)
            {
                i__3 = *n;
                for(i__ = jj + 1; i__ <= i__3; ++i__)
                {
                    i__4 = i__ + (jj - j) * ldwork;
                    i__5 = i__ + jj * a_dim1;
                    work[i__4].real = a[i__5].real;
                    work[i__4].imag = a[i__5].imag; // , expr subst
                    i__4 = i__ + jj * a_dim1;
                    a[i__4].real = 0.f;
                    a[i__4].imag = 0.f; // , expr subst
                    /* L30: */
                }
                /* L40: */
            }
            /* Compute current block column of inv(A). */
            if(j + jb <= *n)
            {
                i__2 = *n - j - jb + 1;
                q__1.real = -1.f;
                q__1.imag = -0.f; // , expr subst
                aocl_blas_cgemm("No transpose", "No transpose", n, &jb, &i__2, &q__1,
                                &a[(j + jb) * a_dim1 + 1], lda, &work[j + jb], &ldwork, &c_b2,
                                &a[j * a_dim1 + 1], lda);
            }
            aocl_blas_ctrsm("Right", "Lower", "No transpose", "Unit", n, &jb, &c_b2, &work[j],
                            &ldwork, &a[j * a_dim1 + 1], lda);
            /* L50: */
        }
    }
    /* Apply column interchanges. */
    for(j = *n - 1; j >= 1; --j)
    {
        jp = ipiv[j];
        if(jp != j)
        {
            aocl_blas_cswap(n, &a[j * a_dim1 + 1], &c__1, &a[jp * a_dim1 + 1], &c__1);
        }
        /* L60: */
    }
    r__1 = aocl_lapack_sroundup_lwork(&iws);
    work[1].real = r__1;
    work[1].imag = 0.f; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGETRI */
}
/* cgetri_ */
