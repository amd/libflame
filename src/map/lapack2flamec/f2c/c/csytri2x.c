/* ../netlib/csytri2x.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {1.f, 0.f};
static scomplex c_b2 = {0.f, 0.f};
/* > \brief \b CSYTRI2X */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CSYTRI2X + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytri2
 * x.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytri2
 * x.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytri2
 * x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N, NB */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ), WORK( N+NB+1,* ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYTRI2X computes the inverse of a real symmetric indefinite matrix */
/* > A using the factorization A = U*D*U**T or A = L*D*L**T computed by */
/* > CSYTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are stored */
/* > as an upper or lower triangular matrix. */
/* > = 'U': Upper triangular, form is A = U*D*U**T;
 */
/* > = 'L': Lower triangular, form is A = L*D*L**T. */
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
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the NNB diagonal matrix D and the multipliers */
/* > used to obtain the factor U or L as computed by CSYTRF. */
/* > */
/* > On exit, if INFO = 0, the (symmetric) inverse of the original */
/* > matrix. If UPLO = 'U', the upper triangular part of the */
/* > inverse is formed and the part of A below the diagonal is not */
/* > referenced;
if UPLO = 'L' the lower triangular part of the */
/* > inverse is formed and the part of A above the diagonal is */
/* > not referenced. */
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
/* > Details of the interchanges and the NNB structure of D */
/* > as determined by CSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (N+NB+1,NB+3) */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > Block size */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, D(i,i) = 0;
the matrix is singular and its */
/* > inverse could not be computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date June 2017 */
/* > \ingroup complexSYcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void csytri2x_(char *uplo, aocl_int_t *n, scomplex *a, aocl_int_t *lda, aocl_int_t *ipiv,
               scomplex *work, aocl_int_t *nb, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_csytri2x(uplo, n, a, lda, ipiv, work, nb, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t nb_64 = *nb;
    aocl_int64_t info_64 = *info;

    aocl_lapack_csytri2x(uplo, &n_64, a, &lda_64, ipiv, work, &nb_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_csytri2x(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                          aocl_int_t *ipiv, scomplex *work, aocl_int64_t *nb, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "csytri2x inputs: uplo %c, n %lld, lda %lld, nb %lld", *uplo, *n, *lda,
             *nb);
#else
    snprintf(buffer, 256, "csytri2x inputs: uplo %c, n %d, lda %d, nb %d", *uplo, *n, *lda, *nb);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, work_dim1, work_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    scomplex q__1, q__2, q__3;
    /* Builtin functions */
    void c_div(scomplex *, scomplex *, scomplex *);
    /* Local variables */
    scomplex d__;
    aocl_int64_t i__, j, k;
    scomplex t, ak;
    aocl_int64_t u11, ip, nnb, cut;
    scomplex akp1;
    aocl_int64_t invd;
    scomplex akkp1;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t iinfo;
    aocl_int64_t count;
    logical upper;
    scomplex u01_i_j__, u11_i_j__;
    scomplex u01_ip1_j__, u11_ip1_j__;
    /* -- LAPACK computational routine (version 3.7.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2017 */
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
    work_dim1 = *n + *nb + 1;
    work_offset = 1 + work_dim1;
    work -= work_offset;
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
    else if(*lda < fla_max(1, *n))
    {
        *info = -4;
    }
    /* Quick return if possible */
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CSYTRI2X", &i__1, (ftnlen)8);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Convert A */
    /* Workspace got Non-diag elements of D */
    aocl_lapack_csyconv(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[work_offset], &iinfo);
    /* Check that the diagonal matrix D is nonsingular. */
    if(upper)
    {
        /* Upper triangular storage: examine D from bottom to top */
        for(*info = *n; *info >= 1; --(*info))
        {
            i__1 = *info + *info * a_dim1;
            if(ipiv[*info] > 0 && (a[i__1].real == 0.f && a[i__1].imag == 0.f))
            {
                AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
                return;
            }
        }
    }
    else
    {
        /* Lower triangular storage: examine D from top to bottom. */
        i__1 = *n;
        for(*info = 1; *info <= i__1; ++(*info))
        {
            i__2 = *info + *info * a_dim1;
            if(ipiv[*info] > 0 && (a[i__2].real == 0.f && a[i__2].imag == 0.f))
            {
                AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
                return;
            }
        }
    }
    *info = 0;
    /* Splitting Workspace */
    /* U01 is a block (N,NB+1) */
    /* The first element of U01 is in WORK(1,1) */
    /* U11 is a block (NB+1,NB+1) */
    /* The first element of U11 is in WORK(N+1,1) */
    u11 = *n;
    /* INVD is a block (N,2) */
    /* The first element of INVD is in WORK(1,INVD) */
    invd = *nb + 2;
    if(upper)
    {
        /* invA = P * inv(U**T)*inv(D)*inv(U)*P**T. */
        aocl_lapack_ctrtri(uplo, "U", n, &a[a_offset], lda, info);
        /* inv(D) and inv(D)*inv(U) */
        k = 1;
        while(k <= *n)
        {
            if(ipiv[k] > 0)
            {
                /* 1 x 1 diagonal NNB */
                i__1 = k + invd * work_dim1;
                c_div(&q__1, &c_b1, &a[k + k * a_dim1]);
                work[i__1].real = q__1.real;
                work[i__1].imag = q__1.imag; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                work[i__1].real = 0.f;
                work[i__1].imag = 0.f; // , expr subst
                ++k;
            }
            else
            {
                /* 2 x 2 diagonal NNB */
                i__1 = k + 1 + work_dim1;
                t.real = work[i__1].real;
                t.imag = work[i__1].imag; // , expr subst
                c_div(&q__1, &a[k + k * a_dim1], &t);
                ak.real = q__1.real;
                ak.imag = q__1.imag; // , expr subst
                c_div(&q__1, &a[k + 1 + (k + 1) * a_dim1], &t);
                akp1.real = q__1.real;
                akp1.imag = q__1.imag; // , expr subst
                c_div(&q__1, &work[k + 1 + work_dim1], &t);
                akkp1.real = q__1.real;
                akkp1.imag = q__1.imag; // , expr subst
                q__3.real = ak.real * akp1.real - ak.imag * akp1.imag;
                q__3.imag = ak.real * akp1.imag + ak.imag * akp1.real; // , expr subst
                q__2.real = q__3.real - 1.f;
                q__2.imag = q__3.imag - 0.f; // , expr subst
                q__1.real = t.real * q__2.real - t.imag * q__2.imag;
                q__1.imag = t.real * q__2.imag + t.imag * q__2.real; // , expr subst
                d__.real = q__1.real;
                d__.imag = q__1.imag; // , expr subst
                i__1 = k + invd * work_dim1;
                c_div(&q__1, &akp1, &d__);
                work[i__1].real = q__1.real;
                work[i__1].imag = q__1.imag; // , expr subst
                i__1 = k + 1 + (invd + 1) * work_dim1;
                c_div(&q__1, &ak, &d__);
                work[i__1].real = q__1.real;
                work[i__1].imag = q__1.imag; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                q__2.real = -akkp1.real;
                q__2.imag = -akkp1.imag; // , expr subst
                c_div(&q__1, &q__2, &d__);
                work[i__1].real = q__1.real;
                work[i__1].imag = q__1.imag; // , expr subst
                i__1 = k + 1 + invd * work_dim1;
                q__2.real = -akkp1.real;
                q__2.imag = -akkp1.imag; // , expr subst
                c_div(&q__1, &q__2, &d__);
                work[i__1].real = q__1.real;
                work[i__1].imag = q__1.imag; // , expr subst
                k += 2;
            }
        }
        /* inv(U**T) = (inv(U))**T */
        /* inv(U**T)*inv(D)*inv(U) */
        cut = *n;
        while(cut > 0)
        {
            nnb = *nb;
            if(cut <= nnb)
            {
                nnb = cut;
            }
            else
            {
                count = 0;
                /* count negative elements, */
                i__1 = cut;
                for(i__ = cut + 1 - nnb; i__ <= i__1; ++i__)
                {
                    if(ipiv[i__] < 0)
                    {
                        ++count;
                    }
                }
                /* need a even number for a clear cut */
                if(count % 2 == 1)
                {
                    ++nnb;
                }
            }
            cut -= nnb;
            /* U01 Block */
            i__1 = cut;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = nnb;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = i__ + j * work_dim1;
                    i__4 = i__ + (cut + j) * a_dim1;
                    work[i__3].real = a[i__4].real;
                    work[i__3].imag = a[i__4].imag; // , expr subst
                }
            }
            /* U11 Block */
            i__1 = nnb;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = u11 + i__ + i__ * work_dim1;
                work[i__2].real = 1.f;
                work[i__2].imag = 0.f; // , expr subst
                i__2 = i__ - 1;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    work[i__3].real = 0.f;
                    work[i__3].imag = 0.f; // , expr subst
                }
                i__2 = nnb;
                for(j = i__ + 1; j <= i__2; ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    i__4 = cut + i__ + (cut + j) * a_dim1;
                    work[i__3].real = a[i__4].real;
                    work[i__3].imag = a[i__4].imag; // , expr subst
                }
            }
            /* invD*U01 */
            i__ = 1;
            while(i__ <= cut)
            {
                if(ipiv[i__] > 0)
                {
                    i__1 = nnb;
                    for(j = 1; j <= i__1; ++j)
                    {
                        i__2 = i__ + j * work_dim1;
                        i__3 = i__ + invd * work_dim1;
                        i__4 = i__ + j * work_dim1;
                        q__1.real = work[i__3].real * work[i__4].real - work[i__3].imag * work[i__4].imag;
                        q__1.imag = work[i__3].real * work[i__4].imag
                                 + work[i__3].imag * work[i__4].real; // , expr subst
                        work[i__2].real = q__1.real;
                        work[i__2].imag = q__1.imag; // , expr subst
                    }
                    ++i__;
                }
                else
                {
                    i__1 = nnb;
                    for(j = 1; j <= i__1; ++j)
                    {
                        i__2 = i__ + j * work_dim1;
                        u01_i_j__.real = work[i__2].real;
                        u01_i_j__.imag = work[i__2].imag; // , expr subst
                        i__2 = i__ + 1 + j * work_dim1;
                        u01_ip1_j__.real = work[i__2].real;
                        u01_ip1_j__.imag = work[i__2].imag; // , expr subst
                        i__2 = i__ + j * work_dim1;
                        i__3 = i__ + invd * work_dim1;
                        q__2.real = work[i__3].real * u01_i_j__.real - work[i__3].imag * u01_i_j__.imag;
                        q__2.imag = work[i__3].real * u01_i_j__.imag
                                 + work[i__3].imag * u01_i_j__.real; // , expr subst
                        i__4 = i__ + (invd + 1) * work_dim1;
                        q__3.real = work[i__4].real * u01_ip1_j__.real - work[i__4].imag * u01_ip1_j__.imag;
                        q__3.imag = work[i__4].real * u01_ip1_j__.imag
                                 + work[i__4].imag * u01_ip1_j__.real; // , expr subst
                        q__1.real = q__2.real + q__3.real;
                        q__1.imag = q__2.imag + q__3.imag; // , expr subst
                        work[i__2].real = q__1.real;
                        work[i__2].imag = q__1.imag; // , expr subst
                        i__2 = i__ + 1 + j * work_dim1;
                        i__3 = i__ + 1 + invd * work_dim1;
                        q__2.real = work[i__3].real * u01_i_j__.real - work[i__3].imag * u01_i_j__.imag;
                        q__2.imag = work[i__3].real * u01_i_j__.imag
                                 + work[i__3].imag * u01_i_j__.real; // , expr subst
                        i__4 = i__ + 1 + (invd + 1) * work_dim1;
                        q__3.real = work[i__4].real * u01_ip1_j__.real - work[i__4].imag * u01_ip1_j__.imag;
                        q__3.imag = work[i__4].real * u01_ip1_j__.imag
                                 + work[i__4].imag * u01_ip1_j__.real; // , expr subst
                        q__1.real = q__2.real + q__3.real;
                        q__1.imag = q__2.imag + q__3.imag; // , expr subst
                        work[i__2].real = q__1.real;
                        work[i__2].imag = q__1.imag; // , expr subst
                    }
                    i__ += 2;
                }
            }
            /* invD1*U11 */
            i__ = 1;
            while(i__ <= nnb)
            {
                if(ipiv[cut + i__] > 0)
                {
                    i__1 = nnb;
                    for(j = i__; j <= i__1; ++j)
                    {
                        i__2 = u11 + i__ + j * work_dim1;
                        i__3 = cut + i__ + invd * work_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        q__1.real = work[i__3].real * work[i__4].real - work[i__3].imag * work[i__4].imag;
                        q__1.imag = work[i__3].real * work[i__4].imag
                                 + work[i__3].imag * work[i__4].real; // , expr subst
                        work[i__2].real = q__1.real;
                        work[i__2].imag = q__1.imag; // , expr subst
                    }
                    ++i__;
                }
                else
                {
                    i__1 = nnb;
                    for(j = i__; j <= i__1; ++j)
                    {
                        i__2 = u11 + i__ + j * work_dim1;
                        u11_i_j__.real = work[i__2].real;
                        u11_i_j__.imag = work[i__2].imag; // , expr subst
                        i__2 = u11 + i__ + 1 + j * work_dim1;
                        u11_ip1_j__.real = work[i__2].real;
                        u11_ip1_j__.imag = work[i__2].imag; // , expr subst
                        i__2 = u11 + i__ + j * work_dim1;
                        i__3 = cut + i__ + invd * work_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        q__2.real = work[i__3].real * work[i__4].real - work[i__3].imag * work[i__4].imag;
                        q__2.imag = work[i__3].real * work[i__4].imag
                                 + work[i__3].imag * work[i__4].real; // , expr subst
                        i__5 = cut + i__ + (invd + 1) * work_dim1;
                        i__6 = u11 + i__ + 1 + j * work_dim1;
                        q__3.real = work[i__5].real * work[i__6].real - work[i__5].imag * work[i__6].imag;
                        q__3.imag = work[i__5].real * work[i__6].imag
                                 + work[i__5].imag * work[i__6].real; // , expr subst
                        q__1.real = q__2.real + q__3.real;
                        q__1.imag = q__2.imag + q__3.imag; // , expr subst
                        work[i__2].real = q__1.real;
                        work[i__2].imag = q__1.imag; // , expr subst
                        i__2 = u11 + i__ + 1 + j * work_dim1;
                        i__3 = cut + i__ + 1 + invd * work_dim1;
                        q__2.real = work[i__3].real * u11_i_j__.real - work[i__3].imag * u11_i_j__.imag;
                        q__2.imag = work[i__3].real * u11_i_j__.imag
                                 + work[i__3].imag * u11_i_j__.real; // , expr subst
                        i__4 = cut + i__ + 1 + (invd + 1) * work_dim1;
                        q__3.real = work[i__4].real * u11_ip1_j__.real - work[i__4].imag * u11_ip1_j__.imag;
                        q__3.imag = work[i__4].real * u11_ip1_j__.imag
                                 + work[i__4].imag * u11_ip1_j__.real; // , expr subst
                        q__1.real = q__2.real + q__3.real;
                        q__1.imag = q__2.imag + q__3.imag; // , expr subst
                        work[i__2].real = q__1.real;
                        work[i__2].imag = q__1.imag; // , expr subst
                    }
                    i__ += 2;
                }
            }
            /* U11**T*invD1*U11->U11 */
            i__1 = *n + *nb + 1;
            aocl_blas_ctrmm("L", "U", "T", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut + 1) * a_dim1],
                            lda, &work[u11 + 1 + work_dim1], &i__1);
            i__1 = nnb;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = nnb;
                for(j = i__; j <= i__2; ++j)
                {
                    i__3 = cut + i__ + (cut + j) * a_dim1;
                    i__4 = u11 + i__ + j * work_dim1;
                    a[i__3].real = work[i__4].real;
                    a[i__3].imag = work[i__4].imag; // , expr subst
                }
            }
            /* U01**T*invD*U01->A(CUT+I,CUT+J) */
            i__1 = *n + *nb + 1;
            i__2 = *n + *nb + 1;
            aocl_blas_cgemm("T", "N", &nnb, &nnb, &cut, &c_b1, &a[(cut + 1) * a_dim1 + 1], lda,
                            &work[work_offset], &i__1, &c_b2, &work[u11 + 1 + work_dim1], &i__2);
            /* U11 = U11**T*invD1*U11 + U01**T*invD*U01 */
            i__1 = nnb;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = nnb;
                for(j = i__; j <= i__2; ++j)
                {
                    i__3 = cut + i__ + (cut + j) * a_dim1;
                    i__4 = cut + i__ + (cut + j) * a_dim1;
                    i__5 = u11 + i__ + j * work_dim1;
                    q__1.real = a[i__4].real + work[i__5].real;
                    q__1.imag = a[i__4].imag + work[i__5].imag; // , expr subst
                    a[i__3].real = q__1.real;
                    a[i__3].imag = q__1.imag; // , expr subst
                }
            }
            /* U01 = U00**T*invD0*U01 */
            i__1 = *n + *nb + 1;
            aocl_blas_ctrmm("L", uplo, "T", "U", &cut, &nnb, &c_b1, &a[a_offset], lda,
                            &work[work_offset], &i__1);
            /* Update U01 */
            i__1 = cut;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = nnb;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = i__ + (cut + j) * a_dim1;
                    i__4 = i__ + j * work_dim1;
                    a[i__3].real = work[i__4].real;
                    a[i__3].imag = work[i__4].imag; // , expr subst
                }
            }
            /* Next Block */
        }
        /* Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T */
        i__ = 1;
        while(i__ <= *n)
        {
            if(ipiv[i__] > 0)
            {
                ip = ipiv[i__];
                if(i__ < ip)
                {
                    aocl_lapack_csyswapr(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if(i__ > ip)
                {
                    aocl_lapack_csyswapr(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
            }
            else
            {
                ip = -ipiv[i__];
                ++i__;
                if(i__ - 1 < ip)
                {
                    i__1 = i__ - 1;
                    aocl_lapack_csyswapr(uplo, n, &a[a_offset], lda, &i__1, &ip);
                }
                if(i__ - 1 > ip)
                {
                    i__1 = i__ - 1;
                    aocl_lapack_csyswapr(uplo, n, &a[a_offset], lda, &ip, &i__1);
                }
            }
            ++i__;
        }
    }
    else
    {
        /* LOWER... */
        /* invA = P * inv(U**T)*inv(D)*inv(U)*P**T. */
        aocl_lapack_ctrtri(uplo, "U", n, &a[a_offset], lda, info);
        /* inv(D) and inv(D)*inv(U) */
        k = *n;
        while(k >= 1)
        {
            if(ipiv[k] > 0)
            {
                /* 1 x 1 diagonal NNB */
                i__1 = k + invd * work_dim1;
                c_div(&q__1, &c_b1, &a[k + k * a_dim1]);
                work[i__1].real = q__1.real;
                work[i__1].imag = q__1.imag; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                work[i__1].real = 0.f;
                work[i__1].imag = 0.f; // , expr subst
                --k;
            }
            else
            {
                /* 2 x 2 diagonal NNB */
                i__1 = k - 1 + work_dim1;
                t.real = work[i__1].real;
                t.imag = work[i__1].imag; // , expr subst
                c_div(&q__1, &a[k - 1 + (k - 1) * a_dim1], &t);
                ak.real = q__1.real;
                ak.imag = q__1.imag; // , expr subst
                c_div(&q__1, &a[k + k * a_dim1], &t);
                akp1.real = q__1.real;
                akp1.imag = q__1.imag; // , expr subst
                c_div(&q__1, &work[k - 1 + work_dim1], &t);
                akkp1.real = q__1.real;
                akkp1.imag = q__1.imag; // , expr subst
                q__3.real = ak.real * akp1.real - ak.imag * akp1.imag;
                q__3.imag = ak.real * akp1.imag + ak.imag * akp1.real; // , expr subst
                q__2.real = q__3.real - 1.f;
                q__2.imag = q__3.imag - 0.f; // , expr subst
                q__1.real = t.real * q__2.real - t.imag * q__2.imag;
                q__1.imag = t.real * q__2.imag + t.imag * q__2.real; // , expr subst
                d__.real = q__1.real;
                d__.imag = q__1.imag; // , expr subst
                i__1 = k - 1 + invd * work_dim1;
                c_div(&q__1, &akp1, &d__);
                work[i__1].real = q__1.real;
                work[i__1].imag = q__1.imag; // , expr subst
                i__1 = k + invd * work_dim1;
                c_div(&q__1, &ak, &d__);
                work[i__1].real = q__1.real;
                work[i__1].imag = q__1.imag; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                q__2.real = -akkp1.real;
                q__2.imag = -akkp1.imag; // , expr subst
                c_div(&q__1, &q__2, &d__);
                work[i__1].real = q__1.real;
                work[i__1].imag = q__1.imag; // , expr subst
                i__1 = k - 1 + (invd + 1) * work_dim1;
                q__2.real = -akkp1.real;
                q__2.imag = -akkp1.imag; // , expr subst
                c_div(&q__1, &q__2, &d__);
                work[i__1].real = q__1.real;
                work[i__1].imag = q__1.imag; // , expr subst
                k += -2;
            }
        }
        /* inv(U**T) = (inv(U))**T */
        /* inv(U**T)*inv(D)*inv(U) */
        cut = 0;
        while(cut < *n)
        {
            nnb = *nb;
            if(cut + nnb >= *n)
            {
                nnb = *n - cut;
            }
            else
            {
                count = 0;
                /* count negative elements, */
                i__1 = cut + nnb;
                for(i__ = cut + 1; i__ <= i__1; ++i__)
                {
                    if(ipiv[i__] < 0)
                    {
                        ++count;
                    }
                }
                /* need a even number for a clear cut */
                if(count % 2 == 1)
                {
                    ++nnb;
                }
            }
            /* L21 Block */
            i__1 = *n - cut - nnb;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = nnb;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = i__ + j * work_dim1;
                    i__4 = cut + nnb + i__ + (cut + j) * a_dim1;
                    work[i__3].real = a[i__4].real;
                    work[i__3].imag = a[i__4].imag; // , expr subst
                }
            }
            /* L11 Block */
            i__1 = nnb;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = u11 + i__ + i__ * work_dim1;
                work[i__2].real = 1.f;
                work[i__2].imag = 0.f; // , expr subst
                i__2 = nnb;
                for(j = i__ + 1; j <= i__2; ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    work[i__3].real = 0.f;
                    work[i__3].imag = 0.f; // , expr subst
                }
                i__2 = i__ - 1;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    i__4 = cut + i__ + (cut + j) * a_dim1;
                    work[i__3].real = a[i__4].real;
                    work[i__3].imag = a[i__4].imag; // , expr subst
                }
            }
            /* invD*L21 */
            i__ = *n - cut - nnb;
            while(i__ >= 1)
            {
                if(ipiv[cut + nnb + i__] > 0)
                {
                    i__1 = nnb;
                    for(j = 1; j <= i__1; ++j)
                    {
                        i__2 = i__ + j * work_dim1;
                        i__3 = cut + nnb + i__ + invd * work_dim1;
                        i__4 = i__ + j * work_dim1;
                        q__1.real = work[i__3].real * work[i__4].real - work[i__3].imag * work[i__4].imag;
                        q__1.imag = work[i__3].real * work[i__4].imag
                                 + work[i__3].imag * work[i__4].real; // , expr subst
                        work[i__2].real = q__1.real;
                        work[i__2].imag = q__1.imag; // , expr subst
                    }
                    --i__;
                }
                else
                {
                    i__1 = nnb;
                    for(j = 1; j <= i__1; ++j)
                    {
                        i__2 = i__ + j * work_dim1;
                        u01_i_j__.real = work[i__2].real;
                        u01_i_j__.imag = work[i__2].imag; // , expr subst
                        i__2 = i__ - 1 + j * work_dim1;
                        u01_ip1_j__.real = work[i__2].real;
                        u01_ip1_j__.imag = work[i__2].imag; // , expr subst
                        i__2 = i__ + j * work_dim1;
                        i__3 = cut + nnb + i__ + invd * work_dim1;
                        q__2.real = work[i__3].real * u01_i_j__.real - work[i__3].imag * u01_i_j__.imag;
                        q__2.imag = work[i__3].real * u01_i_j__.imag
                                 + work[i__3].imag * u01_i_j__.real; // , expr subst
                        i__4 = cut + nnb + i__ + (invd + 1) * work_dim1;
                        q__3.real = work[i__4].real * u01_ip1_j__.real - work[i__4].imag * u01_ip1_j__.imag;
                        q__3.imag = work[i__4].real * u01_ip1_j__.imag
                                 + work[i__4].imag * u01_ip1_j__.real; // , expr subst
                        q__1.real = q__2.real + q__3.real;
                        q__1.imag = q__2.imag + q__3.imag; // , expr subst
                        work[i__2].real = q__1.real;
                        work[i__2].imag = q__1.imag; // , expr subst
                        i__2 = i__ - 1 + j * work_dim1;
                        i__3 = cut + nnb + i__ - 1 + (invd + 1) * work_dim1;
                        q__2.real = work[i__3].real * u01_i_j__.real - work[i__3].imag * u01_i_j__.imag;
                        q__2.imag = work[i__3].real * u01_i_j__.imag
                                 + work[i__3].imag * u01_i_j__.real; // , expr subst
                        i__4 = cut + nnb + i__ - 1 + invd * work_dim1;
                        q__3.real = work[i__4].real * u01_ip1_j__.real - work[i__4].imag * u01_ip1_j__.imag;
                        q__3.imag = work[i__4].real * u01_ip1_j__.imag
                                 + work[i__4].imag * u01_ip1_j__.real; // , expr subst
                        q__1.real = q__2.real + q__3.real;
                        q__1.imag = q__2.imag + q__3.imag; // , expr subst
                        work[i__2].real = q__1.real;
                        work[i__2].imag = q__1.imag; // , expr subst
                    }
                    i__ += -2;
                }
            }
            /* invD1*L11 */
            i__ = nnb;
            while(i__ >= 1)
            {
                if(ipiv[cut + i__] > 0)
                {
                    i__1 = nnb;
                    for(j = 1; j <= i__1; ++j)
                    {
                        i__2 = u11 + i__ + j * work_dim1;
                        i__3 = cut + i__ + invd * work_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        q__1.real = work[i__3].real * work[i__4].real - work[i__3].imag * work[i__4].imag;
                        q__1.imag = work[i__3].real * work[i__4].imag
                                 + work[i__3].imag * work[i__4].real; // , expr subst
                        work[i__2].real = q__1.real;
                        work[i__2].imag = q__1.imag; // , expr subst
                    }
                    --i__;
                }
                else
                {
                    i__1 = nnb;
                    for(j = 1; j <= i__1; ++j)
                    {
                        i__2 = u11 + i__ + j * work_dim1;
                        u11_i_j__.real = work[i__2].real;
                        u11_i_j__.imag = work[i__2].imag; // , expr subst
                        i__2 = u11 + i__ - 1 + j * work_dim1;
                        u11_ip1_j__.real = work[i__2].real;
                        u11_ip1_j__.imag = work[i__2].imag; // , expr subst
                        i__2 = u11 + i__ + j * work_dim1;
                        i__3 = cut + i__ + invd * work_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        q__2.real = work[i__3].real * work[i__4].real - work[i__3].imag * work[i__4].imag;
                        q__2.imag = work[i__3].real * work[i__4].imag
                                 + work[i__3].imag * work[i__4].real; // , expr subst
                        i__5 = cut + i__ + (invd + 1) * work_dim1;
                        q__3.real = work[i__5].real * u11_ip1_j__.real - work[i__5].imag * u11_ip1_j__.imag;
                        q__3.imag = work[i__5].real * u11_ip1_j__.imag
                                 + work[i__5].imag * u11_ip1_j__.real; // , expr subst
                        q__1.real = q__2.real + q__3.real;
                        q__1.imag = q__2.imag + q__3.imag; // , expr subst
                        work[i__2].real = q__1.real;
                        work[i__2].imag = q__1.imag; // , expr subst
                        i__2 = u11 + i__ - 1 + j * work_dim1;
                        i__3 = cut + i__ - 1 + (invd + 1) * work_dim1;
                        q__2.real = work[i__3].real * u11_i_j__.real - work[i__3].imag * u11_i_j__.imag;
                        q__2.imag = work[i__3].real * u11_i_j__.imag
                                 + work[i__3].imag * u11_i_j__.real; // , expr subst
                        i__4 = cut + i__ - 1 + invd * work_dim1;
                        q__3.real = work[i__4].real * u11_ip1_j__.real - work[i__4].imag * u11_ip1_j__.imag;
                        q__3.imag = work[i__4].real * u11_ip1_j__.imag
                                 + work[i__4].imag * u11_ip1_j__.real; // , expr subst
                        q__1.real = q__2.real + q__3.real;
                        q__1.imag = q__2.imag + q__3.imag; // , expr subst
                        work[i__2].real = q__1.real;
                        work[i__2].imag = q__1.imag; // , expr subst
                    }
                    i__ += -2;
                }
            }
            /* L11**T*invD1*L11->L11 */
            i__1 = *n + *nb + 1;
            aocl_blas_ctrmm("L", uplo, "T", "U", &nnb, &nnb, &c_b1,
                            &a[cut + 1 + (cut + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1],
                            &i__1);
            i__1 = nnb;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = i__;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = cut + i__ + (cut + j) * a_dim1;
                    i__4 = u11 + i__ + j * work_dim1;
                    a[i__3].real = work[i__4].real;
                    a[i__3].imag = work[i__4].imag; // , expr subst
                }
            }
            if(cut + nnb < *n)
            {
                /* L21**T*invD2*L21->A(CUT+I,CUT+J) */
                i__1 = *n - nnb - cut;
                i__2 = *n + *nb + 1;
                i__3 = *n + *nb + 1;
                aocl_blas_cgemm("T", "N", &nnb, &nnb, &i__1, &c_b1,
                                &a[cut + nnb + 1 + (cut + 1) * a_dim1], lda, &work[work_offset],
                                &i__2, &c_b2, &work[u11 + 1 + work_dim1], &i__3);
                /* L11 = L11**T*invD1*L11 + U01**T*invD*U01 */
                i__1 = nnb;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = i__;
                    for(j = 1; j <= i__2; ++j)
                    {
                        i__3 = cut + i__ + (cut + j) * a_dim1;
                        i__4 = cut + i__ + (cut + j) * a_dim1;
                        i__5 = u11 + i__ + j * work_dim1;
                        q__1.real = a[i__4].real + work[i__5].real;
                        q__1.imag = a[i__4].imag + work[i__5].imag; // , expr subst
                        a[i__3].real = q__1.real;
                        a[i__3].imag = q__1.imag; // , expr subst
                    }
                }
                /* L01 = L22**T*invD2*L21 */
                i__1 = *n - nnb - cut;
                i__2 = *n + *nb + 1;
                aocl_blas_ctrmm("L", uplo, "T", "U", &i__1, &nnb, &c_b1,
                                &a[cut + nnb + 1 + (cut + nnb + 1) * a_dim1], lda,
                                &work[work_offset], &i__2);
                /* Update L21 */
                i__1 = *n - cut - nnb;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = nnb;
                    for(j = 1; j <= i__2; ++j)
                    {
                        i__3 = cut + nnb + i__ + (cut + j) * a_dim1;
                        i__4 = i__ + j * work_dim1;
                        a[i__3].real = work[i__4].real;
                        a[i__3].imag = work[i__4].imag; // , expr subst
                    }
                }
            }
            else
            {
                /* L11 = L11**T*invD1*L11 */
                i__1 = nnb;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = i__;
                    for(j = 1; j <= i__2; ++j)
                    {
                        i__3 = cut + i__ + (cut + j) * a_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        a[i__3].real = work[i__4].real;
                        a[i__3].imag = work[i__4].imag; // , expr subst
                    }
                }
            }
            /* Next Block */
            cut += nnb;
        }
        /* Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T */
        i__ = *n;
        while(i__ >= 1)
        {
            if(ipiv[i__] > 0)
            {
                ip = ipiv[i__];
                if(i__ < ip)
                {
                    aocl_lapack_csyswapr(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if(i__ > ip)
                {
                    aocl_lapack_csyswapr(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
            }
            else
            {
                ip = -ipiv[i__];
                if(i__ < ip)
                {
                    aocl_lapack_csyswapr(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if(i__ > ip)
                {
                    aocl_lapack_csyswapr(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
                --i__;
            }
            --i__;
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CSYTRI2X */
}
/* csytri2x_ */
