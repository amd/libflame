/* ../netlib/v3.9.0/chb2st_kernels.f -- translated by f2c (version 20160102). You must link the
 resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or
 Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place,
 with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief \b CHB2ST_KERNELS */
/* @generated from zhb2st_kernels.f, fortran z -> c, Wed Dec 7 08:22:40 2016 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHB2ST_KERNELS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chb2st_
 * kernels.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chb2st_
 * kernels.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chb2st_
 * kernels.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHB2ST_KERNELS( UPLO, WANTZ, TTYPE, */
/* ST, ED, SWEEP, N, NB, IB, */
/* A, LDA, V, TAU, LDVT, WORK) */
/* IMPLICIT NONE */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* LOGICAL WANTZ */
/* INTEGER TTYPE, ST, ED, SWEEP, N, NB, IB, LDA, LDVT */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), V( * ), */
/* TAU( * ), WORK( * ) */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHB2ST_KERNELS is an internal routine used by the CHETRD_HB2ST */
/* > subroutine. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* > WANTZ is LOGICAL which indicate if Eigenvalue are requested or both */
/* > Eigenvalue/Eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] TTYPE */
/* > \verbatim */
/* > TTYPE is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] ST */
/* > \verbatim */
/* > ST is INTEGER */
/* > internal parameter for indices. */
/* > \endverbatim */
/* > */
/* > \param[in] ED */
/* > \verbatim */
/* > ED is INTEGER */
/* > internal parameter for indices. */
/* > \endverbatim */
/* > */
/* > \param[in] SWEEP */
/* > \verbatim */
/* > SWEEP is INTEGER */
/* > internal parameter for indices. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER. The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER. The size of the band. */
/* > \endverbatim */
/* > */
/* > \param[in] IB */
/* > \verbatim */
/* > IB is INTEGER. */
/* > \endverbatim */
/* > */
/* > \param[in, out] A */
/* > \verbatim */
/* > A is COMPLEX array. A pointer to the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER. The leading dimension of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* > V is COMPLEX array, dimension 2*n if eigenvalues only are */
/* > requested or to be queried for vectors. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX array, dimension (2*n). */
/* > The scalar factors of the Householder reflectors are stored */
/* > in this array. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* > LDVT is INTEGER. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array. Workspace of size nb. */
/* > \endverbatim */
/* > */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Implemented by Azzam Haidar. */
/* > */
/* > All details are available on technical report, SC11, SC13 papers. */
/* > */
/* > Azzam Haidar, Hatem Ltaief, and Jack Dongarra. */
/* > Parallel reduction to condensed forms for symmetric eigenvalue problems */
/* > using aggregated fine-grained and memory-aware kernels. In Proceedings */
/* > of 2011 International Conference for High Performance Computing, */
/* > Networking, Storage and Analysis (SC '11), New York, NY, USA, */
/* > Article 8 , 11 pages. */
/* > http://doi.acm.org/10.1145/2063384.2063394 */
/* > */
/* > A. Haidar, J. Kurzak, P. Luszczek, 2013. */
/* > An improved parallel singular value algorithm and its implementation */
/* > for multicore hardware, In Proceedings of 2013 International Conference */
/* > for High Performance Computing, Networking, Storage and Analysis (SC '13). */
/* > Denver, Colorado, USA, 2013. */
/* > Article 90, 12 pages. */
/* > http://doi.acm.org/10.1145/2503210.2503292 */
/* > */
/* > A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra. */
/* > A novel hybrid CPU-GPU generalized eigensolver for electronic structure */
/* > calculations based on fine-grained memory aware tasks. */
/* > International Journal of High Performance Computing Applications. */
/* > Volume 28 Issue 2, Pages 196-209, May 2014. */
/* > http://hpc.sagepub.com/content/28/2/196 */
/* > */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void chb2st_kernels_(char *uplo, logical *wantz, aocl_int_t *ttype, aocl_int_t *st, aocl_int_t *ed,
                     aocl_int_t *sweep, aocl_int_t *n, aocl_int_t *nb, aocl_int_t *ib, scomplex *a,
                     aocl_int_t *lda, scomplex *v, scomplex *tau, aocl_int_t *ldvt, scomplex *work)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_chb2st_kernels(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt,
                               work);
#else
    aocl_int64_t ttype_64 = *ttype;
    aocl_int64_t st_64 = *st;
    aocl_int64_t ed_64 = *ed;
    aocl_int64_t sweep_64 = *sweep;
    aocl_int64_t n_64 = *n;
    aocl_int64_t nb_64 = *nb;
    aocl_int64_t ib_64 = *ib;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldvt_64 = *ldvt;

    aocl_lapack_chb2st_kernels(uplo, wantz, &ttype_64, &st_64, &ed_64, &sweep_64, &n_64, &nb_64,
                               &ib_64, a, &lda_64, v, tau, &ldvt_64, work);
#endif
}

void aocl_lapack_chb2st_kernels(char *uplo, logical *wantz, aocl_int64_t *ttype, aocl_int64_t *st,
                                aocl_int64_t *ed, aocl_int64_t *sweep, aocl_int64_t *n,
                                aocl_int64_t *nb, aocl_int64_t *ib, scomplex *a, aocl_int64_t *lda,
                                scomplex *v, scomplex *tau, aocl_int64_t *ldvt, scomplex *work)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "chb2st_kernels inputs: uplo %c, ttype %lld, st %lld, ed %lld, sweep %lld, n %lld, nb "
             "%lld, ib %lld, lda %lld, ldvt %lld",
             *uplo, *ttype, *st, *ed, *sweep, *n, *nb, *ib, *lda, *ldvt);
#else
    snprintf(buffer, 256,
             "chb2st_kernels inputs: uplo %c, ttype %d, st %d, ed %d, sweep %d, n %d, nb %d, ib "
             "%d, lda %d, ldvt %d",
             *uplo, *ttype, *st, *ed, *sweep, *n, *nb, *ib, *lda, *ldvt);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3;
    scomplex q__1;
    /* Builtin functions */
    void r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    aocl_int64_t i__, j1, j2, lm, ln;
    scomplex ctmp;
    aocl_int64_t dpos, vpos;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    logical upper;
    aocl_int64_t ofdpos, taupos;
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. External Functions .. */
    /* .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --v;
    --tau;
    --work;
    /* Function Body */
    upper = lsame_(uplo, "U", 1, 1);
    if(upper)
    {
        dpos = (*nb << 1) + 1;
        ofdpos = *nb << 1;
    }
    else
    {
        dpos = 1;
        ofdpos = 2;
    }
    /* Upper case */
    if(upper)
    {
        if(*wantz)
        {
            vpos = (*sweep - 1) % 2 * *n + *st;
            taupos = (*sweep - 1) % 2 * *n + *st;
        }
        else
        {
            vpos = (*sweep - 1) % 2 * *n + *st;
            taupos = (*sweep - 1) % 2 * *n + *st;
        }
        if(*ttype == 1)
        {
            lm = *ed - *st + 1;
            i__1 = vpos;
            v[i__1].r = 1.f;
            v[i__1].i = 0.f; // , expr subst
            i__1 = lm - 1;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = vpos + i__;
                r_cnjg(&q__1, &a[ofdpos - i__ + (*st + i__) * a_dim1]);
                v[i__2].r = q__1.r;
                v[i__2].i = q__1.i; // , expr subst
                i__2 = ofdpos - i__ + (*st + i__) * a_dim1;
                a[i__2].r = 0.f;
                a[i__2].i = 0.f; // , expr subst
                /* L10: */
            }
            r_cnjg(&q__1, &a[ofdpos + *st * a_dim1]);
            ctmp.r = q__1.r;
            ctmp.i = q__1.i; // , expr subst
            aocl_lapack_clarfg(&lm, &ctmp, &v[vpos + 1], &c__1, &tau[taupos]);
            i__1 = ofdpos + *st * a_dim1;
            a[i__1].r = ctmp.r;
            a[i__1].i = ctmp.i; // , expr subst
            lm = *ed - *st + 1;
            r_cnjg(&q__1, &tau[taupos]);
            i__1 = *lda - 1;
            aocl_lapack_clarfy(uplo, &lm, &v[vpos], &c__1, &q__1, &a[dpos + *st * a_dim1], &i__1,
                               &work[1]);
        }
        if(*ttype == 3)
        {
            lm = *ed - *st + 1;
            r_cnjg(&q__1, &tau[taupos]);
            i__1 = *lda - 1;
            aocl_lapack_clarfy(uplo, &lm, &v[vpos], &c__1, &q__1, &a[dpos + *st * a_dim1], &i__1,
                               &work[1]);
        }
        if(*ttype == 2)
        {
            j1 = *ed + 1;
            /* Computing MIN */
            i__1 = *ed + *nb;
            j2 = fla_min(i__1, *n);
            ln = *ed - *st + 1;
            lm = j2 - j1 + 1;
            if(lm > 0)
            {
                r_cnjg(&q__1, &tau[taupos]);
                i__1 = *lda - 1;
                aocl_lapack_clarfx("Left", &ln, &lm, &v[vpos], &q__1, &a[dpos - *nb + j1 * a_dim1],
                                   &i__1, &work[1]);
                if(*wantz)
                {
                    vpos = (*sweep - 1) % 2 * *n + j1;
                    taupos = (*sweep - 1) % 2 * *n + j1;
                }
                else
                {
                    vpos = (*sweep - 1) % 2 * *n + j1;
                    taupos = (*sweep - 1) % 2 * *n + j1;
                }
                i__1 = vpos;
                v[i__1].r = 1.f;
                v[i__1].i = 0.f; // , expr subst
                i__1 = lm - 1;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = vpos + i__;
                    r_cnjg(&q__1, &a[dpos - *nb - i__ + (j1 + i__) * a_dim1]);
                    v[i__2].r = q__1.r;
                    v[i__2].i = q__1.i; // , expr subst
                    i__2 = dpos - *nb - i__ + (j1 + i__) * a_dim1;
                    a[i__2].r = 0.f;
                    a[i__2].i = 0.f; // , expr subst
                    /* L30: */
                }
                r_cnjg(&q__1, &a[dpos - *nb + j1 * a_dim1]);
                ctmp.r = q__1.r;
                ctmp.i = q__1.i; // , expr subst
                aocl_lapack_clarfg(&lm, &ctmp, &v[vpos + 1], &c__1, &tau[taupos]);
                i__1 = dpos - *nb + j1 * a_dim1;
                a[i__1].r = ctmp.r;
                a[i__1].i = ctmp.i; // , expr subst
                i__1 = ln - 1;
                i__2 = *lda - 1;
                aocl_lapack_clarfx("Right", &i__1, &lm, &v[vpos], &tau[taupos],
                                   &a[dpos - *nb + 1 + j1 * a_dim1], &i__2, &work[1]);
            }
        }
        /* Lower case */
    }
    else
    {
        if(*wantz)
        {
            vpos = (*sweep - 1) % 2 * *n + *st;
            taupos = (*sweep - 1) % 2 * *n + *st;
        }
        else
        {
            vpos = (*sweep - 1) % 2 * *n + *st;
            taupos = (*sweep - 1) % 2 * *n + *st;
        }
        if(*ttype == 1)
        {
            lm = *ed - *st + 1;
            i__1 = vpos;
            v[i__1].r = 1.f;
            v[i__1].i = 0.f; // , expr subst
            i__1 = lm - 1;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = vpos + i__;
                i__3 = ofdpos + i__ + (*st - 1) * a_dim1;
                v[i__2].r = a[i__3].r;
                v[i__2].i = a[i__3].i; // , expr subst
                i__2 = ofdpos + i__ + (*st - 1) * a_dim1;
                a[i__2].r = 0.f;
                a[i__2].i = 0.f; // , expr subst
                /* L20: */
            }
            aocl_lapack_clarfg(&lm, &a[ofdpos + (*st - 1) * a_dim1], &v[vpos + 1], &c__1,
                               &tau[taupos]);
            lm = *ed - *st + 1;
            r_cnjg(&q__1, &tau[taupos]);
            i__1 = *lda - 1;
            aocl_lapack_clarfy(uplo, &lm, &v[vpos], &c__1, &q__1, &a[dpos + *st * a_dim1], &i__1,
                               &work[1]);
        }
        if(*ttype == 3)
        {
            lm = *ed - *st + 1;
            r_cnjg(&q__1, &tau[taupos]);
            i__1 = *lda - 1;
            aocl_lapack_clarfy(uplo, &lm, &v[vpos], &c__1, &q__1, &a[dpos + *st * a_dim1], &i__1,
                               &work[1]);
        }
        if(*ttype == 2)
        {
            j1 = *ed + 1;
            /* Computing MIN */
            i__1 = *ed + *nb;
            j2 = fla_min(i__1, *n);
            ln = *ed - *st + 1;
            lm = j2 - j1 + 1;
            if(lm > 0)
            {
                i__1 = *lda - 1;
                aocl_lapack_clarfx("Right", &lm, &ln, &v[vpos], &tau[taupos],
                                   &a[dpos + *nb + *st * a_dim1], &i__1, &work[1]);
                if(*wantz)
                {
                    vpos = (*sweep - 1) % 2 * *n + j1;
                    taupos = (*sweep - 1) % 2 * *n + j1;
                }
                else
                {
                    vpos = (*sweep - 1) % 2 * *n + j1;
                    taupos = (*sweep - 1) % 2 * *n + j1;
                }
                i__1 = vpos;
                v[i__1].r = 1.f;
                v[i__1].i = 0.f; // , expr subst
                i__1 = lm - 1;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = vpos + i__;
                    i__3 = dpos + *nb + i__ + *st * a_dim1;
                    v[i__2].r = a[i__3].r;
                    v[i__2].i = a[i__3].i; // , expr subst
                    i__2 = dpos + *nb + i__ + *st * a_dim1;
                    a[i__2].r = 0.f;
                    a[i__2].i = 0.f; // , expr subst
                    /* L40: */
                }
                aocl_lapack_clarfg(&lm, &a[dpos + *nb + *st * a_dim1], &v[vpos + 1], &c__1,
                                   &tau[taupos]);
                i__1 = ln - 1;
                r_cnjg(&q__1, &tau[taupos]);
                i__2 = *lda - 1;
                aocl_lapack_clarfx("Left", &lm, &i__1, &v[vpos], &q__1,
                                   &a[dpos + *nb - 1 + (*st + 1) * a_dim1], &i__2, &work[1]);
            }
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* END OF CHB2ST_KERNELS */
}
/* chb2st_kernels__ */
