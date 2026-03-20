/*
    Copyright (c) 2019-2023 Advanced Micro Devices, Inc.
*/
/* ../netlib/zgeqp3.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
static aocl_int64_t c__3 = 3;
static aocl_int64_t c__2 = 2;
/* > \brief \b ZGEQP3 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGEQP3 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeqp3.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeqp3.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeqp3.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, RWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER JPVT( * ) */
/* DOUBLE PRECISION RWORK( * ) */
/* COMPLEX*16 A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEQP3 computes a QR factorization with column pivoting of a */
/* > matrix A: A*P = Q*R using Level 3 BLAS. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, the upper triangle of the array contains the */
/* > fla_min(M,N)-by-N upper trapezoidal matrix R;
the elements below */
/* > the diagonal, together with the array TAU, represent the */
/* > unitary matrix Q as a product of fla_min(M,N) elementary */
/* > reflectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] JPVT */
/* > \verbatim */
/* > JPVT is INTEGER array, dimension (N) */
/* > On entry, if JPVT(J).ne.0, the J-th column of A is permuted */
/* > to the front of A*P (a leading column);
if JPVT(J)=0, */
/* > the J-th column of A is a free column. */
/* > On exit, if JPVT(J)=K, then the J-th column of A*P was the */
/* > the K-th column of A. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX*16 array, dimension (fla_min(M,N)) */
/* > The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO=0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= N+1. */
/* > For optimal performance LWORK >= ( N+1 )*NB, where NB */
/* > is the optimal blocksize. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is DOUBLE PRECISION array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup complex16GEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix Q is represented as a product of elementary reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(k), where k = fla_min(m,n). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**H */
/* > */
/* > where tau is a scomplex scalar, and v is a real/scomplex vector */
/* > with v(1:i-1) = 0 and v(i) = 1;
v(i+1:m) is stored on exit in */
/* > A(i+1:m,i), and tau in TAU(i). */
/* > \endverbatim */
/* > \par Contributors: */
/* ================== */
/* > */
/* > G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain */
/* > X. Sun, Computer Science Dept., Duke University, USA */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zgeqp3_(aocl_int_t *m, aocl_int_t *n, dcomplex *a, aocl_int_t *lda, aocl_int_t *jpvt,
             dcomplex *tau, dcomplex *work, aocl_int_t *lwork, doublereal *rwork,
             aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zgeqp3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zgeqp3(&m_64, &n_64, a, &lda_64, jpvt, tau, work, &lwork_64, rwork, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_zgeqp3(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                        aocl_int_t *jpvt, dcomplex *tau, dcomplex *work,
                        aocl_int64_t *lwork, doublereal *rwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgeqp3 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", lwork %" FLA_IS
                      "",
                      *m, *n, *lda, *lwork);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3;
    dcomplex z__1;
    /* Local variables */
    aocl_int64_t j, jb, na, nb, sm, sn, nx, fjb, iws, nfxd, nbmin, minmn, minws;
    aocl_int64_t topbmn, sminmn;
    aocl_int64_t lwkopt;
    logical lquery;
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test input arguments */
    /* ==================== */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --jpvt;
    --tau;
    --work;
    --rwork;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -4;
    }
    if(*info == 0)
    {
        minmn = fla_min(*m, *n);
        if(minmn == 0)
        {
            iws = 1;
            lwkopt = 1;
        }
        else
        {
            iws = *n + 1;
#ifdef FLA_ENABLE_AMD_OPT
            nb = FLA_ZGEQP3_BLOCK_SMALL_THRESH;
#else
            nb = aocl_lapack_ilaenv(&c__1, "ZGEQRF", " ", m, n, &c_n1, &c_n1);
#endif
            lwkopt = (*n + 1) * nb;
        }
        z__1.real = (doublereal)lwkopt;
        z__1.imag = 0.; // , expr subst
        work[1].real = z__1.real;
        work[1].imag = z__1.imag; // , expr subst
        if(*lwork < iws && !lquery)
        {
            *info = -8;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("ZGEQP3", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Move initial columns up front. */
    nfxd = 1;
    i__1 = *n;
    for(j = 1; j <= i__1; ++j)
    {
        if(jpvt[j] != 0)
        {
            if(j != nfxd)
            {
                aocl_blas_zswap(m, &a[j * a_dim1 + 1], &c__1, &a[nfxd * a_dim1 + 1], &c__1);
                jpvt[j] = jpvt[nfxd];
                jpvt[nfxd] = (aocl_int_t)(j);
            }
            else
            {
                jpvt[j] = (aocl_int_t)(j);
            }
            ++nfxd;
        }
        else
        {
            jpvt[j] = (aocl_int_t)(j);
        }
        /* L10: */
    }
    --nfxd;
    /* Factorize fixed columns */
    /* ======================= */
    /* Compute the QR factorization of fixed columns and update */
    /* remaining columns. */
    if(nfxd > 0)
    {
        na = fla_min(*m, nfxd);
        /* CC CALL ZGEQR2( M, NA, A, LDA, TAU, WORK, INFO ) */
        aocl_lapack_zgeqrf(m, &na, &a[a_offset], lda, &tau[1], &work[1], lwork, info);
        /* Computing MAX */
        i__1 = iws;
        i__2 = (integer)work[1].real; // , expr subst
        iws = fla_max(i__1, i__2);
        if(na < *n)
        {
            /* CC CALL ZUNM2R( 'Left', 'Conjugate Transpose', M, N-NA, */
            /* CC $ NA, A, LDA, TAU, A( 1, NA+1 ), LDA, WORK, */
            /* CC $ INFO ) */
            i__1 = *n - na;
            aocl_lapack_zunmqr("Left", "Conjugate Transpose", m, &i__1, &na, &a[a_offset], lda,
                               &tau[1], &a[(na + 1) * a_dim1 + 1], lda, &work[1], lwork, info);
            /* Computing MAX */
            i__1 = iws;
            i__2 = (integer)work[1].real; // , expr subst
            iws = fla_max(i__1, i__2);
        }
    }
    /* Factorize free columns */
    /* ====================== */
    if(nfxd < minmn)
    {
        sm = *m - nfxd;
        sn = *n - nfxd;
        sminmn = minmn - nfxd;
        /* Determine the block size. */
#ifdef FLA_ENABLE_AMD_OPT
        nb = FLA_ZGEQP3_BLOCK_SMALL_THRESH;
#else
        nb = aocl_lapack_ilaenv(&c__1, "ZGEQRF", " ", &sm, &sn, &c_n1, &c_n1);
#endif
        nbmin = 2;
        nx = 0;
        if(nb > 1 && nb < sminmn)
        {
            /* Determine when to cross over from blocked to unblocked code. */
            /* Computing MAX */
            i__1 = 0;
            i__2 = aocl_lapack_ilaenv(&c__3, "ZGEQRF", " ", &sm, &sn, &c_n1, &c_n1); // , expr subst
            nx = fla_max(i__1, i__2);
            if(nx < sminmn)
            {
                /* Determine if workspace is large enough for blocked code. */
                minws = (sn + 1) * nb;
                iws = fla_max(iws, minws);
                if(*lwork < minws)
                {
                    /* Not enough workspace to use optimal NB: Reduce NB and */
                    /* determine the minimum value of NB. */
                    nb = *lwork / (sn + 1);
                    /* Computing MAX */
                    i__1 = 2;
                    i__2 = aocl_lapack_ilaenv(&c__2, "ZGEQRF", " ", &sm, &sn, &c_n1,
                                              &c_n1); // , expr subst
                    nbmin = fla_max(i__1, i__2);
                }
            }
        }
        /* Initialize partial column norms. The first N elements of work */
        /* store the exact column norms. */
        i__1 = *n;
        for(j = nfxd + 1; j <= i__1; ++j)
        {
            rwork[j] = aocl_blas_dznrm2(&sm, &a[nfxd + 1 + j * a_dim1], &c__1);
            rwork[*n + j] = rwork[j];
            /* L20: */
        }
        if(nb >= nbmin && nb < sminmn && nx < sminmn)
        {
            /* Use blocked code initially. */
            j = nfxd + 1;
            /* Compute factorization: while loop. */
            topbmn = minmn - nx;
        L30:
            if(j <= topbmn)
            {
                /* Computing MIN */
                i__1 = nb;
                i__2 = topbmn - j + 1; // , expr subst
                jb = fla_min(i__1, i__2);
                /* Factorize JB columns among columns J:N. */
                i__1 = *n - j + 1;
                i__2 = j - 1;
                i__3 = *n - j + 1;
                aocl_lapack_zlaqps(m, &i__1, &i__2, &jb, &fjb, &a[j * a_dim1 + 1], lda, &jpvt[j],
                                   &tau[j], &rwork[j], &rwork[*n + j], &work[1], &work[jb + 1],
                                   &i__3);
                j += fjb;
                goto L30;
            }
        }
        else
        {
            j = nfxd + 1;
        }
        /* Use unblocked code to factor the last or only block. */
        if(j <= minmn)
        {
            i__1 = *n - j + 1;
            i__2 = j - 1;
            aocl_lapack_zlaqp2(m, &i__1, &i__2, &a[j * a_dim1 + 1], lda, &jpvt[j], &tau[j],
                               &rwork[j], &rwork[*n + j], &work[1]);
        }
    }
    z__1.real = (doublereal)lwkopt;
    z__1.imag = 0.; // , expr subst
    work[1].real = z__1.real;
    work[1].imag = z__1.imag; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZGEQP3 */
}
/* zgeqp3_ */
