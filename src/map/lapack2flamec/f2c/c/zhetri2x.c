/* ../netlib/zhetri2x.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static dcomplex c_b1 = {1., 0.};
static dcomplex c_b2 = {0., 0.};
/* > \brief \b ZHETRI2X */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHETRI2X + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetri2
 * x.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetri2
 * x.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetri2
 * x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHETRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N, NB */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), WORK( N+NB+1,* ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHETRI2X computes the inverse of a COMPLEX*16 Hermitian indefinite matrix */
/* > A using the factorization A = U*D*U**H or A = L*D*L**H computed by */
/* > ZHETRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are stored */
/* > as an upper or lower triangular matrix. */
/* > = 'U': Upper triangular, form is A = U*D*U**H;
 */
/* > = 'L': Lower triangular, form is A = L*D*L**H. */
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
/* > On entry, the NNB diagonal matrix D and the multipliers */
/* > used to obtain the factor U or L as computed by ZHETRF. */
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
/* > as determined by ZHETRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (N+NB+1,NB+3) */
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
/* > \ingroup complex16HEcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zhetri2x_(char *uplo, aocl_int_t *n, dcomplex *a, aocl_int_t *lda, aocl_int_t *ipiv,
               dcomplex *work, aocl_int_t *nb, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zhetri2x(uplo, n, a, lda, ipiv, work, nb, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t nb_64 = *nb;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zhetri2x(uplo, &n_64, a, &lda_64, ipiv, work, &nb_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_zhetri2x(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                          aocl_int_t *ipiv, dcomplex *work, aocl_int64_t *nb,
                          aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhetri2x inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS ", nb %" FLA_IS "",
                      *uplo, *n, *lda, *nb);

    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, work_dim1, work_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    dcomplex z__1, z__2, z__3;
    /* Builtin functions */
    double z_abs(dcomplex *);
    void z_div(dcomplex *, dcomplex *, dcomplex *),
        d_cnjg(dcomplex *, dcomplex *);
    /* Local variables */
    dcomplex d__;
    aocl_int64_t i__, j, k;
    dcomplex t, ak;
    aocl_int64_t u11, ip, nnb, cut;
    dcomplex akp1;
    aocl_int64_t invd;
    dcomplex akkp1;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t iinfo;
    aocl_int64_t count;
    logical upper;
    dcomplex u01_i_j__, u11_i_j__;
    dcomplex u01_ip1_j__, u11_ip1_j__;
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
        aocl_blas_xerbla("ZHETRI2X", &i__1, (ftnlen)8);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Convert A */
    /* Workspace got Non-diag elements of D */
    aocl_lapack_zsyconv(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[work_offset], &iinfo);
    /* Check that the diagonal matrix D is nonsingular. */
    if(upper)
    {
        /* Upper triangular storage: examine D from bottom to top */
        for(*info = *n; *info >= 1; --(*info))
        {
            i__1 = *info + *info * a_dim1;
            if(ipiv[*info] > 0 && (a[i__1].real == 0. && a[i__1].imag == 0.))
            {
                AOCL_DTL_TRACE_LOG_EXIT
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
            if(ipiv[*info] > 0 && (a[i__2].real == 0. && a[i__2].imag == 0.))
            {
                AOCL_DTL_TRACE_LOG_EXIT
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
        /* invA = P * inv(U**H)*inv(D)*inv(U)*P**H. */
        aocl_lapack_ztrtri(uplo, "U", n, &a[a_offset], lda, info);
        /* inv(D) and inv(D)*inv(U) */
        k = 1;
        while(k <= *n)
        {
            if(ipiv[k] > 0)
            {
                /* 1 x 1 diagonal NNB */
                i__1 = k + invd * work_dim1;
                i__2 = k + k * a_dim1;
                d__1 = 1. / a[i__2].real;
                work[i__1].real = d__1;
                work[i__1].imag = 0.; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                work[i__1].real = 0.;
                work[i__1].imag = 0.; // , expr subst
                ++k;
            }
            else
            {
                /* 2 x 2 diagonal NNB */
                d__1 = z_abs(&work[k + 1 + work_dim1]);
                t.real = d__1;
                t.imag = 0.; // , expr subst
                i__1 = k + k * a_dim1;
                d__1 = a[i__1].real;
                z__2.real = d__1;
                z__2.imag = 0.; // , expr subst
                z_div(&z__1, &z__2, &t);
                ak.real = z__1.real;
                ak.imag = z__1.imag; // , expr subst
                i__1 = k + 1 + (k + 1) * a_dim1;
                d__1 = a[i__1].real;
                z__2.real = d__1;
                z__2.imag = 0.; // , expr subst
                z_div(&z__1, &z__2, &t);
                akp1.real = z__1.real;
                akp1.imag = z__1.imag; // , expr subst
                z_div(&z__1, &work[k + 1 + work_dim1], &t);
                akkp1.real = z__1.real;
                akkp1.imag = z__1.imag; // , expr subst
                z__3.real = ak.real * akp1.real - ak.imag * akp1.imag;
                z__3.imag = ak.real * akp1.imag + ak.imag * akp1.real; // , expr subst
                z__2.real = z__3.real - 1.;
                z__2.imag = z__3.imag; // , expr subst
                z__1.real = t.real * z__2.real - t.imag * z__2.imag;
                z__1.imag = t.real * z__2.imag + t.imag * z__2.real; // , expr subst
                d__.real = z__1.real;
                d__.imag = z__1.imag; // , expr subst
                i__1 = k + invd * work_dim1;
                z_div(&z__1, &akp1, &d__);
                work[i__1].real = z__1.real;
                work[i__1].imag = z__1.imag; // , expr subst
                i__1 = k + 1 + (invd + 1) * work_dim1;
                z_div(&z__1, &ak, &d__);
                work[i__1].real = z__1.real;
                work[i__1].imag = z__1.imag; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                z__2.real = -akkp1.real;
                z__2.imag = -akkp1.imag; // , expr subst
                z_div(&z__1, &z__2, &d__);
                work[i__1].real = z__1.real;
                work[i__1].imag = z__1.imag; // , expr subst
                i__1 = k + 1 + invd * work_dim1;
                d_cnjg(&z__1, &work[k + (invd + 1) * work_dim1]);
                work[i__1].real = z__1.real;
                work[i__1].imag = z__1.imag; // , expr subst
                k += 2;
            }
        }
        /* inv(U**H) = (inv(U))**H */
        /* inv(U**H)*inv(D)*inv(U) */
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
                work[i__2].real = 1.;
                work[i__2].imag = 0.; // , expr subst
                i__2 = i__ - 1;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    work[i__3].real = 0.;
                    work[i__3].imag = 0.; // , expr subst
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
                        z__1.real = work[i__3].real * work[i__4].real - work[i__3].imag * work[i__4].imag;
                        z__1.imag = work[i__3].real * work[i__4].imag
                                 + work[i__3].imag * work[i__4].real; // , expr subst
                        work[i__2].real = z__1.real;
                        work[i__2].imag = z__1.imag; // , expr subst
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
                        z__2.real = work[i__3].real * u01_i_j__.real - work[i__3].imag * u01_i_j__.imag;
                        z__2.imag = work[i__3].real * u01_i_j__.imag
                                 + work[i__3].imag * u01_i_j__.real; // , expr subst
                        i__4 = i__ + (invd + 1) * work_dim1;
                        z__3.real = work[i__4].real * u01_ip1_j__.real - work[i__4].imag * u01_ip1_j__.imag;
                        z__3.imag = work[i__4].real * u01_ip1_j__.imag
                                 + work[i__4].imag * u01_ip1_j__.real; // , expr subst
                        z__1.real = z__2.real + z__3.real;
                        z__1.imag = z__2.imag + z__3.imag; // , expr subst
                        work[i__2].real = z__1.real;
                        work[i__2].imag = z__1.imag; // , expr subst
                        i__2 = i__ + 1 + j * work_dim1;
                        i__3 = i__ + 1 + invd * work_dim1;
                        z__2.real = work[i__3].real * u01_i_j__.real - work[i__3].imag * u01_i_j__.imag;
                        z__2.imag = work[i__3].real * u01_i_j__.imag
                                 + work[i__3].imag * u01_i_j__.real; // , expr subst
                        i__4 = i__ + 1 + (invd + 1) * work_dim1;
                        z__3.real = work[i__4].real * u01_ip1_j__.real - work[i__4].imag * u01_ip1_j__.imag;
                        z__3.imag = work[i__4].real * u01_ip1_j__.imag
                                 + work[i__4].imag * u01_ip1_j__.real; // , expr subst
                        z__1.real = z__2.real + z__3.real;
                        z__1.imag = z__2.imag + z__3.imag; // , expr subst
                        work[i__2].real = z__1.real;
                        work[i__2].imag = z__1.imag; // , expr subst
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
                        z__1.real = work[i__3].real * work[i__4].real - work[i__3].imag * work[i__4].imag;
                        z__1.imag = work[i__3].real * work[i__4].imag
                                 + work[i__3].imag * work[i__4].real; // , expr subst
                        work[i__2].real = z__1.real;
                        work[i__2].imag = z__1.imag; // , expr subst
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
                        z__2.real = work[i__3].real * work[i__4].real - work[i__3].imag * work[i__4].imag;
                        z__2.imag = work[i__3].real * work[i__4].imag
                                 + work[i__3].imag * work[i__4].real; // , expr subst
                        i__5 = cut + i__ + (invd + 1) * work_dim1;
                        i__6 = u11 + i__ + 1 + j * work_dim1;
                        z__3.real = work[i__5].real * work[i__6].real - work[i__5].imag * work[i__6].imag;
                        z__3.imag = work[i__5].real * work[i__6].imag
                                 + work[i__5].imag * work[i__6].real; // , expr subst
                        z__1.real = z__2.real + z__3.real;
                        z__1.imag = z__2.imag + z__3.imag; // , expr subst
                        work[i__2].real = z__1.real;
                        work[i__2].imag = z__1.imag; // , expr subst
                        i__2 = u11 + i__ + 1 + j * work_dim1;
                        i__3 = cut + i__ + 1 + invd * work_dim1;
                        z__2.real = work[i__3].real * u11_i_j__.real - work[i__3].imag * u11_i_j__.imag;
                        z__2.imag = work[i__3].real * u11_i_j__.imag
                                 + work[i__3].imag * u11_i_j__.real; // , expr subst
                        i__4 = cut + i__ + 1 + (invd + 1) * work_dim1;
                        z__3.real = work[i__4].real * u11_ip1_j__.real - work[i__4].imag * u11_ip1_j__.imag;
                        z__3.imag = work[i__4].real * u11_ip1_j__.imag
                                 + work[i__4].imag * u11_ip1_j__.real; // , expr subst
                        z__1.real = z__2.real + z__3.real;
                        z__1.imag = z__2.imag + z__3.imag; // , expr subst
                        work[i__2].real = z__1.real;
                        work[i__2].imag = z__1.imag; // , expr subst
                    }
                    i__ += 2;
                }
            }
            /* U11**H*invD1*U11->U11 */
            i__1 = *n + *nb + 1;
            aocl_blas_ztrmm("L", "U", "C", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut + 1) * a_dim1],
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
            /* U01**H*invD*U01->A(CUT+I,CUT+J) */
            i__1 = *n + *nb + 1;
            i__2 = *n + *nb + 1;
            aocl_blas_zgemm("C", "N", &nnb, &nnb, &cut, &c_b1, &a[(cut + 1) * a_dim1 + 1], lda,
                            &work[work_offset], &i__1, &c_b2, &work[u11 + 1 + work_dim1], &i__2);
            /* U11 = U11**H*invD1*U11 + U01**H*invD*U01 */
            i__1 = nnb;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = nnb;
                for(j = i__; j <= i__2; ++j)
                {
                    i__3 = cut + i__ + (cut + j) * a_dim1;
                    i__4 = cut + i__ + (cut + j) * a_dim1;
                    i__5 = u11 + i__ + j * work_dim1;
                    z__1.real = a[i__4].real + work[i__5].real;
                    z__1.imag = a[i__4].imag + work[i__5].imag; // , expr subst
                    a[i__3].real = z__1.real;
                    a[i__3].imag = z__1.imag; // , expr subst
                }
            }
            /* U01 = U00**H*invD0*U01 */
            i__1 = *n + *nb + 1;
            aocl_blas_ztrmm("L", uplo, "C", "U", &cut, &nnb, &c_b1, &a[a_offset], lda,
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
        /* Apply PERMUTATIONS P and P**H: P * inv(U**H)*inv(D)*inv(U) *P**H */
        i__ = 1;
        while(i__ <= *n)
        {
            if(ipiv[i__] > 0)
            {
                ip = ipiv[i__];
                if(i__ < ip)
                {
                    aocl_lapack_zheswapr(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if(i__ > ip)
                {
                    aocl_lapack_zheswapr(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
            }
            else
            {
                ip = -ipiv[i__];
                ++i__;
                if(i__ - 1 < ip)
                {
                    i__1 = i__ - 1;
                    aocl_lapack_zheswapr(uplo, n, &a[a_offset], lda, &i__1, &ip);
                }
                if(i__ - 1 > ip)
                {
                    i__1 = i__ - 1;
                    aocl_lapack_zheswapr(uplo, n, &a[a_offset], lda, &ip, &i__1);
                }
            }
            ++i__;
        }
    }
    else
    {
        /* LOWER... */
        /* invA = P * inv(U**H)*inv(D)*inv(U)*P**H. */
        aocl_lapack_ztrtri(uplo, "U", n, &a[a_offset], lda, info);
        /* inv(D) and inv(D)*inv(U) */
        k = *n;
        while(k >= 1)
        {
            if(ipiv[k] > 0)
            {
                /* 1 x 1 diagonal NNB */
                i__1 = k + invd * work_dim1;
                i__2 = k + k * a_dim1;
                d__1 = 1. / a[i__2].real;
                work[i__1].real = d__1;
                work[i__1].imag = 0.; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                work[i__1].real = 0.;
                work[i__1].imag = 0.; // , expr subst
                --k;
            }
            else
            {
                /* 2 x 2 diagonal NNB */
                d__1 = z_abs(&work[k - 1 + work_dim1]);
                t.real = d__1;
                t.imag = 0.; // , expr subst
                i__1 = k - 1 + (k - 1) * a_dim1;
                d__1 = a[i__1].real;
                z__2.real = d__1;
                z__2.imag = 0.; // , expr subst
                z_div(&z__1, &z__2, &t);
                ak.real = z__1.real;
                ak.imag = z__1.imag; // , expr subst
                i__1 = k + k * a_dim1;
                d__1 = a[i__1].real;
                z__2.real = d__1;
                z__2.imag = 0.; // , expr subst
                z_div(&z__1, &z__2, &t);
                akp1.real = z__1.real;
                akp1.imag = z__1.imag; // , expr subst
                z_div(&z__1, &work[k - 1 + work_dim1], &t);
                akkp1.real = z__1.real;
                akkp1.imag = z__1.imag; // , expr subst
                z__3.real = ak.real * akp1.real - ak.imag * akp1.imag;
                z__3.imag = ak.real * akp1.imag + ak.imag * akp1.real; // , expr subst
                z__2.real = z__3.real - 1.;
                z__2.imag = z__3.imag; // , expr subst
                z__1.real = t.real * z__2.real - t.imag * z__2.imag;
                z__1.imag = t.real * z__2.imag + t.imag * z__2.real; // , expr subst
                d__.real = z__1.real;
                d__.imag = z__1.imag; // , expr subst
                i__1 = k - 1 + invd * work_dim1;
                z_div(&z__1, &akp1, &d__);
                work[i__1].real = z__1.real;
                work[i__1].imag = z__1.imag; // , expr subst
                i__1 = k + invd * work_dim1;
                z_div(&z__1, &ak, &d__);
                work[i__1].real = z__1.real;
                work[i__1].imag = z__1.imag; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                z__2.real = -akkp1.real;
                z__2.imag = -akkp1.imag; // , expr subst
                z_div(&z__1, &z__2, &d__);
                work[i__1].real = z__1.real;
                work[i__1].imag = z__1.imag; // , expr subst
                i__1 = k - 1 + (invd + 1) * work_dim1;
                d_cnjg(&z__1, &work[k + (invd + 1) * work_dim1]);
                work[i__1].real = z__1.real;
                work[i__1].imag = z__1.imag; // , expr subst
                k += -2;
            }
        }
        /* inv(U**H) = (inv(U))**H */
        /* inv(U**H)*inv(D)*inv(U) */
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
                work[i__2].real = 1.;
                work[i__2].imag = 0.; // , expr subst
                i__2 = nnb;
                for(j = i__ + 1; j <= i__2; ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    work[i__3].real = 0.;
                    work[i__3].imag = 0.; // , expr subst
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
                        z__1.real = work[i__3].real * work[i__4].real - work[i__3].imag * work[i__4].imag;
                        z__1.imag = work[i__3].real * work[i__4].imag
                                 + work[i__3].imag * work[i__4].real; // , expr subst
                        work[i__2].real = z__1.real;
                        work[i__2].imag = z__1.imag; // , expr subst
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
                        z__2.real = work[i__3].real * u01_i_j__.real - work[i__3].imag * u01_i_j__.imag;
                        z__2.imag = work[i__3].real * u01_i_j__.imag
                                 + work[i__3].imag * u01_i_j__.real; // , expr subst
                        i__4 = cut + nnb + i__ + (invd + 1) * work_dim1;
                        z__3.real = work[i__4].real * u01_ip1_j__.real - work[i__4].imag * u01_ip1_j__.imag;
                        z__3.imag = work[i__4].real * u01_ip1_j__.imag
                                 + work[i__4].imag * u01_ip1_j__.real; // , expr subst
                        z__1.real = z__2.real + z__3.real;
                        z__1.imag = z__2.imag + z__3.imag; // , expr subst
                        work[i__2].real = z__1.real;
                        work[i__2].imag = z__1.imag; // , expr subst
                        i__2 = i__ - 1 + j * work_dim1;
                        i__3 = cut + nnb + i__ - 1 + (invd + 1) * work_dim1;
                        z__2.real = work[i__3].real * u01_i_j__.real - work[i__3].imag * u01_i_j__.imag;
                        z__2.imag = work[i__3].real * u01_i_j__.imag
                                 + work[i__3].imag * u01_i_j__.real; // , expr subst
                        i__4 = cut + nnb + i__ - 1 + invd * work_dim1;
                        z__3.real = work[i__4].real * u01_ip1_j__.real - work[i__4].imag * u01_ip1_j__.imag;
                        z__3.imag = work[i__4].real * u01_ip1_j__.imag
                                 + work[i__4].imag * u01_ip1_j__.real; // , expr subst
                        z__1.real = z__2.real + z__3.real;
                        z__1.imag = z__2.imag + z__3.imag; // , expr subst
                        work[i__2].real = z__1.real;
                        work[i__2].imag = z__1.imag; // , expr subst
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
                        z__1.real = work[i__3].real * work[i__4].real - work[i__3].imag * work[i__4].imag;
                        z__1.imag = work[i__3].real * work[i__4].imag
                                 + work[i__3].imag * work[i__4].real; // , expr subst
                        work[i__2].real = z__1.real;
                        work[i__2].imag = z__1.imag; // , expr subst
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
                        z__2.real = work[i__3].real * work[i__4].real - work[i__3].imag * work[i__4].imag;
                        z__2.imag = work[i__3].real * work[i__4].imag
                                 + work[i__3].imag * work[i__4].real; // , expr subst
                        i__5 = cut + i__ + (invd + 1) * work_dim1;
                        z__3.real = work[i__5].real * u11_ip1_j__.real - work[i__5].imag * u11_ip1_j__.imag;
                        z__3.imag = work[i__5].real * u11_ip1_j__.imag
                                 + work[i__5].imag * u11_ip1_j__.real; // , expr subst
                        z__1.real = z__2.real + z__3.real;
                        z__1.imag = z__2.imag + z__3.imag; // , expr subst
                        work[i__2].real = z__1.real;
                        work[i__2].imag = z__1.imag; // , expr subst
                        i__2 = u11 + i__ - 1 + j * work_dim1;
                        i__3 = cut + i__ - 1 + (invd + 1) * work_dim1;
                        z__2.real = work[i__3].real * u11_i_j__.real - work[i__3].imag * u11_i_j__.imag;
                        z__2.imag = work[i__3].real * u11_i_j__.imag
                                 + work[i__3].imag * u11_i_j__.real; // , expr subst
                        i__4 = cut + i__ - 1 + invd * work_dim1;
                        z__3.real = work[i__4].real * u11_ip1_j__.real - work[i__4].imag * u11_ip1_j__.imag;
                        z__3.imag = work[i__4].real * u11_ip1_j__.imag
                                 + work[i__4].imag * u11_ip1_j__.real; // , expr subst
                        z__1.real = z__2.real + z__3.real;
                        z__1.imag = z__2.imag + z__3.imag; // , expr subst
                        work[i__2].real = z__1.real;
                        work[i__2].imag = z__1.imag; // , expr subst
                    }
                    i__ += -2;
                }
            }
            /* L11**H*invD1*L11->L11 */
            i__1 = *n + *nb + 1;
            aocl_blas_ztrmm("L", uplo, "C", "U", &nnb, &nnb, &c_b1,
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
                /* L21**H*invD2*L21->A(CUT+I,CUT+J) */
                i__1 = *n - nnb - cut;
                i__2 = *n + *nb + 1;
                i__3 = *n + *nb + 1;
                aocl_blas_zgemm("C", "N", &nnb, &nnb, &i__1, &c_b1,
                                &a[cut + nnb + 1 + (cut + 1) * a_dim1], lda, &work[work_offset],
                                &i__2, &c_b2, &work[u11 + 1 + work_dim1], &i__3);
                /* L11 = L11**H*invD1*L11 + U01**H*invD*U01 */
                i__1 = nnb;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = i__;
                    for(j = 1; j <= i__2; ++j)
                    {
                        i__3 = cut + i__ + (cut + j) * a_dim1;
                        i__4 = cut + i__ + (cut + j) * a_dim1;
                        i__5 = u11 + i__ + j * work_dim1;
                        z__1.real = a[i__4].real + work[i__5].real;
                        z__1.imag = a[i__4].imag + work[i__5].imag; // , expr subst
                        a[i__3].real = z__1.real;
                        a[i__3].imag = z__1.imag; // , expr subst
                    }
                }
                /* L01 = L22**H*invD2*L21 */
                i__1 = *n - nnb - cut;
                i__2 = *n + *nb + 1;
                aocl_blas_ztrmm("L", uplo, "C", "U", &i__1, &nnb, &c_b1,
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
                /* L11 = L11**H*invD1*L11 */
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
        /* Apply PERMUTATIONS P and P**H: P * inv(U**H)*inv(D)*inv(U) *P**H */
        i__ = *n;
        while(i__ >= 1)
        {
            if(ipiv[i__] > 0)
            {
                ip = ipiv[i__];
                if(i__ < ip)
                {
                    aocl_lapack_zheswapr(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if(i__ > ip)
                {
                    aocl_lapack_zheswapr(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
            }
            else
            {
                ip = -ipiv[i__];
                if(i__ < ip)
                {
                    aocl_lapack_zheswapr(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if(i__ > ip)
                {
                    aocl_lapack_zheswapr(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
                --i__;
            }
            --i__;
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZHETRI2X */
}
/* zhetri2x_ */
