/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#ifndef FLA_GEQRF_SMALL_CMPLX_AVX2_KERNEL_H
#define FLA_GEQRF_SMALL_CMPLX_AVX2_KERNEL_H

/*! @file fla_geqrf_small_cmplx_avx2_kernel.h
 *  @brief QR kernels for small sizes (complex c/z).
 *  */

#if FLA_ENABLE_AMD_OPT

#define FLA_GEQRF_COMPLEX_SMALL_CAT_(a, b) a##b
#define FLA_GEQRF_COMPLEX_SMALL_CAT(a, b) FLA_GEQRF_COMPLEX_SMALL_CAT_(a, b)

#define c_VW 4
#define c_V __m256
#define c_VLOADU(p) _mm256_loadu_ps((const real *)(p))
#define c_VSTOREU(p, v) _mm256_storeu_ps((real *)(p), v)
#define c_VSET1(x) _mm256_set1_ps(x)
#define c_VZERO() _mm256_setzero_ps()
#define c_VMUL(a, b) _mm256_mul_ps(a, b)
#define c_VADD(a, b) _mm256_add_ps(a, b)
#define c_VSUB(a, b) _mm256_sub_ps(a, b)
#define c_VADDSUB(a, b) _mm256_addsub_ps(a, b)
#define c_VSWAP(v) _mm256_shuffle_ps(v, v, 0xB1)
#define c_VDUPRE(v) _mm256_shuffle_ps(v, v, 0xA0)
#define c_VDUPIM(v) _mm256_shuffle_ps(v, v, 0xF5)
#define c_VHSUM(acc, RE, IM)                             \
    {                                                    \
        __m128 lo_ = _mm256_castps256_ps128(acc);        \
        __m128 hi_ = _mm256_extractf128_ps(acc, 0x1);    \
        __m128 s_ = _mm_add_ps(lo_, hi_);                \
        s_ = _mm_add_ps(s_, _mm_movehl_ps(s_, s_));      \
        RE = _mm_cvtss_f32(s_);                          \
        IM = _mm_cvtss_f32(_mm_shuffle_ps(s_, s_, 0x1)); \
    }

#define z_VW 2
#define z_V __m256d
#define z_VLOADU(p) _mm256_loadu_pd((const doublereal *)(p))
#define z_VSTOREU(p, v) _mm256_storeu_pd((doublereal *)(p), v)
#define z_VSET1(x) _mm256_set1_pd(x)
#define z_VZERO() _mm256_setzero_pd()
#define z_VMUL(a, b) _mm256_mul_pd(a, b)
#define z_VADD(a, b) _mm256_add_pd(a, b)
#define z_VSUB(a, b) _mm256_sub_pd(a, b)
#define z_VADDSUB(a, b) _mm256_addsub_pd(a, b)
#define z_VSWAP(v) _mm256_permute_pd(v, 0x5)
#define z_VDUPRE(v) _mm256_permute_pd(v, 0x0)
#define z_VDUPIM(v) _mm256_permute_pd(v, 0xF)
#define z_VHSUM(acc, RE, IM)                           \
    {                                                  \
        __m128d lo_ = _mm256_castpd256_pd128(acc);     \
        __m128d hi_ = _mm256_extractf128_pd(acc, 0x1); \
        lo_ = _mm_add_pd(lo_, hi_);                    \
        RE = _mm_cvtsd_f64(lo_);                       \
        IM = _mm_cvtsd_f64(_mm_unpackhi_pd(lo_, lo_)); \
    }

/* norm of complex tail (scnrm2/dznrm2 scaled accumulation) */
#define FLA_GEQRF_CMPLX_NRM2(RTYPE, RABS, RSQRT, RES, CXPTR, LEN) \
    {                                                             \
        RTYPE scale_ = (RTYPE)0, ssq_ = (RTYPE)1;                 \
        aocl_int64_t kk_;                                         \
        for(kk_ = 1; kk_ <= (LEN); kk_++)                         \
        {                                                         \
            RTYPE parts_[2];                                      \
            aocl_int64_t pp_;                                     \
            parts_[0] = (CXPTR)[kk_].real;                        \
            parts_[1] = (CXPTR)[kk_].imag;                        \
            for(pp_ = 0; pp_ < 2; pp_++)                          \
            {                                                     \
                RTYPE av_ = RABS(parts_[pp_]);                    \
                if(av_ != (RTYPE)0)                               \
                {                                                 \
                    if(scale_ < av_)                              \
                    {                                             \
                        RTYPE t_ = scale_ / av_;                  \
                        ssq_ = (RTYPE)1 + ssq_ * t_ * t_;         \
                        scale_ = av_;                             \
                    }                                             \
                    else                                          \
                    {                                             \
                        RTYPE t_ = av_ / scale_;                  \
                        ssq_ += t_ * t_;                          \
                    }                                             \
                }                                                 \
            }                                                     \
        }                                                         \
        RES = scale_ * RSQRT(ssq_);                               \
    }

/* beta = -sign(lapy3(alphr, alphi, xnorm), alphr) */
#define FLA_GEQRF_CMPLX_BETA(RTYPE, RABS, RSQRT, BETA, ALPHR, ALPHI, XNORM) \
    {                                                                       \
        RTYPE xa_ = RABS(ALPHR), ya_ = RABS(ALPHI), za_ = RABS(XNORM);      \
        RTYPE w_ = xa_ > ya_ ? xa_ : ya_;                                   \
        RTYPE lp3_;                                                         \
        w_ = w_ > za_ ? w_ : za_;                                           \
        if(w_ == (RTYPE)0)                                                  \
            lp3_ = xa_ + ya_ + za_;                                         \
        else                                                                \
        {                                                                   \
            RTYPE rx_ = xa_ / w_, ry_ = ya_ / w_, rz_ = za_ / w_;           \
            lp3_ = w_ * RSQRT(rx_ * rx_ + ry_ * ry_ + rz_ * rz_);           \
        }                                                                   \
        BETA = (ALPHR >= (RTYPE)0) ? -lp3_ : lp3_;                          \
    }

#define FLA_DEFINE_GEQRF_CMPLX_SMALL(FN, P, CTYPE, RTYPE, LAMCH, RABS, RSQRT)                      \
    int FN(aocl_int64_t *m, aocl_int64_t *n, CTYPE *a, aocl_int64_t *lda, CTYPE *tau, CTYPE *work) \
    {                                                                                              \
        aocl_int64_t i, k, jc, min_m_n, arows, tail, knt;                                          \
        RTYPE xnorm, beta, alphr, alphi, safmin, rsafmn;                                           \
        (void)work;                                                                                \
        /* Adjust pointers */                                                                      \
        a -= (1 + *lda);                                                                           \
        tau--;                                                                                     \
        safmin = LAMCH("S") / LAMCH("E");                                                          \
        rsafmn = (RTYPE)1 / safmin;                                                                \
        min_m_n = fla_min(*m, *n);                                                                 \
        for(i = 1; i <= min_m_n; i++)                                                              \
        {                                                                                          \
            CTYPE *col = &a[i + i * *lda]; /* &A(i,i) */                                           \
            arows = *m - i + 1;                                                                    \
            tail = arows - 1;                                                                      \
            alphr = col[0].real;                                                                   \
            alphi = col[0].imag;                                                                   \
            /* calculate norm of sub-diagonal elements in current column */                        \
            FLA_GEQRF_CMPLX_NRM2(RTYPE, RABS, RSQRT, xnorm, col, tail);                            \
            /* Compute Householder Reflector parameters */                                         \
            if(xnorm == (RTYPE)0 && alphi == (RTYPE)0)                                             \
            {                                                                                      \
                tau[i].real = (RTYPE)0;                                                            \
                tau[i].imag = (RTYPE)0;                                                            \
                continue;                                                                          \
            }                                                                                      \
            FLA_GEQRF_CMPLX_BETA(RTYPE, RABS, RSQRT, beta, alphr, alphi, xnorm);                   \
            /* Scale-up the inputs for small norm */                                               \
            knt = 0;                                                                               \
            if(RABS(beta) < safmin)                                                                \
            {                                                                                      \
                do                                                                                 \
                {                                                                                  \
                    knt++;                                                                         \
                    for(k = 1; k <= tail; k++)                                                     \
                    {                                                                              \
                        col[k].real *= rsafmn;                                                     \
                        col[k].imag *= rsafmn;                                                     \
                    }                                                                              \
                    beta *= rsafmn;                                                                \
                    alphi *= rsafmn;                                                               \
                    alphr *= rsafmn;                                                               \
                } while(RABS(beta) < safmin && knt < 20);                                          \
                FLA_GEQRF_CMPLX_NRM2(RTYPE, RABS, RSQRT, xnorm, col, tail);                        \
                FLA_GEQRF_CMPLX_BETA(RTYPE, RABS, RSQRT, beta, alphr, alphi, xnorm);               \
            }                                                                                      \
            /* Part 1: Compute Householder vector 'v' and tau */                                   \
            tau[i].real = (beta - alphr) / beta;                                                   \
            tau[i].imag = -alphi / beta;                                                           \
            /* scal = 1 / (alpha - beta) */                                                        \
            {                                                                                      \
                RTYPE dr = alphr - beta, di = alphi, sr, si, e_, f_;                               \
                if(RABS(dr) >= RABS(di))                                                           \
                {                                                                                  \
                    e_ = di / dr;                                                                  \
                    f_ = dr + di * e_;                                                             \
                    sr = (RTYPE)1 / f_;                                                            \
                    si = -e_ / f_;                                                                 \
                }                                                                                  \
                else                                                                               \
                {                                                                                  \
                    e_ = dr / di;                                                                  \
                    f_ = dr * e_ + di;                                                             \
                    sr = e_ / f_;                                                                  \
                    si = (RTYPE)(-1) / f_;                                                         \
                }                                                                                  \
                {                                                                                  \
                    FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                             \
                    vsr = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VSET1)(sr);                              \
                    FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                             \
                    vsi = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VSET1)(si);                              \
                    k = 1;                                                                         \
                    for(; k + (FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VW) - 1) <= tail;                   \
                        k += FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VW))                                  \
                    {                                                                              \
                        FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                         \
                        vv = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VLOADU)(&col[k]);                     \
                        FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                         \
                        p1 = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VMUL)(vv, vsr);                       \
                        FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                         \
                        sw = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VSWAP)(vv);                           \
                        FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                         \
                        p2 = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VMUL)(sw, vsi);                       \
                        FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VSTOREU)                                   \
                        (&col[k], FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VADDSUB)(p1, p2));               \
                    }                                                                              \
                    for(; k <= tail; k++)                                                          \
                    {                                                                              \
                        RTYPE vr = col[k].real, vi = col[k].imag;                                  \
                        col[k].real = vr * sr - vi * si;                                           \
                        col[k].imag = vr * si + vi * sr;                                           \
                    }                                                                              \
                }                                                                                  \
            }                                                                                      \
            /* Scale-down beta */                                                                  \
            for(k = 1; k <= knt; k++)                                                              \
                beta *= safmin;                                                                    \
            col[0].real = beta;                                                                    \
            col[0].imag = (RTYPE)0;                                                                \
            /* Part 2: Apply the Householder rotation              */                              \
            /* on the rest of the matrix                           */                              \
            /*    A = A - conj(tau) * v * v**H * A               */                                \
            if(i < *n)                                                                             \
            {                                                                                      \
                RTYPE tcr = tau[i].real, tci = -tau[i].imag;                                       \
                for(jc = i + 1; jc <= *n; jc++)                                                    \
                {                                                                                  \
                    CTYPE *cj = &a[i + jc * *lda]; /* &A(i,jc) */                                  \
                    RTYPE wr, wi, tr, ti;                                                          \
                    /* Compute w = v**H * c_A */                                                   \
                    {                                                                              \
                        FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                         \
                        acc = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VZERO)();                            \
                        k = 1;                                                                     \
                        for(; k + (FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VW) - 1) <= tail;               \
                            k += FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VW))                              \
                        {                                                                          \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                     \
                            cvec = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VLOADU)(&cj[k]);                \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                     \
                            vvec = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VLOADU)(&col[k]);               \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                     \
                            vr_dup = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VDUPRE)(vvec);                \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                     \
                            vi_dup = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VDUPIM)(vvec);                \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                     \
                            c_sw = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VSWAP)(cvec);                   \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                     \
                            aa = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VMUL)(vr_dup, cvec);              \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                     \
                            nvi = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VSUB)(                           \
                                FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VZERO)(), vi_dup);                 \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                     \
                            bb = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VMUL)(nvi, c_sw);                 \
                            acc = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VADD)(                           \
                                acc, FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VADDSUB)(aa, bb));            \
                        }                                                                          \
                        FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VHSUM)(acc, wr, wi);                       \
                        wr += cj[0].real;                                                          \
                        wi += cj[0].imag;                                                          \
                        for(; k <= tail; k++)                                                      \
                        {                                                                          \
                            RTYPE vr = col[k].real, vi = col[k].imag;                              \
                            RTYPE cr = cj[k].real, ci = cj[k].imag;                                \
                            wr += vr * cr + vi * ci;                                               \
                            wi += vr * ci - vi * cr;                                               \
                        }                                                                          \
                    }                                                                              \
                    /* t = conj(tau) * w */                                                        \
                    tr = tcr * wr - tci * wi;                                                      \
                    ti = tcr * wi + tci * wr;                                                      \
                    /* c_A = c_A - v * t */                                                        \
                    cj[0].real -= tr;                                                              \
                    cj[0].imag -= ti;                                                              \
                    {                                                                              \
                        FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                         \
                        tr_b = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VSET1)(tr);                         \
                        FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                         \
                        ti_b = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VSET1)(ti);                         \
                        k = 1;                                                                     \
                        for(; k + (FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VW) - 1) <= tail;               \
                            k += FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VW))                              \
                        {                                                                          \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                     \
                            vvec = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VLOADU)(&col[k]);               \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                     \
                            cvec = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VLOADU)(&cj[k]);                \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                     \
                            aa = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VMUL)(vvec, tr_b);                \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                     \
                            v_sw = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VSWAP)(vvec);                   \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                     \
                            bb = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VMUL)(v_sw, ti_b);                \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _V)                                     \
                            vt = FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VADDSUB)(aa, bb);                 \
                            FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VSTOREU)                               \
                            (&cj[k], FLA_GEQRF_COMPLEX_SMALL_CAT(P, _VSUB)(cvec, vt));             \
                        }                                                                          \
                        for(; k <= tail; k++)                                                      \
                        {                                                                          \
                            RTYPE vr = col[k].real, vi = col[k].imag;                              \
                            cj[k].real -= vr * tr - vi * ti;                                       \
                            cj[k].imag -= vr * ti + vi * tr;                                       \
                        }                                                                          \
                    }                                                                              \
                }                                                                                  \
            }                                                                                      \
        }                                                                                          \
        return 0;                                                                                  \
    }

#endif /* FLA_ENABLE_AMD_OPT */
#endif /* FLA_GEQRF_SMALL_CMPLX_AVX2_KERNEL_H */
