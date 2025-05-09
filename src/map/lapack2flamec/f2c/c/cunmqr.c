/* ../netlib/cunmqr.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
/*
 *  Copyright (c) 2020-2023 Advanced Micro Devices, Inc.  All rights reserved.
 */
#include "FLA_f2c.h" /* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;
/* > \brief \b CUNMQR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CUNMQR + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmqr.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmqr.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmqr.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
/* WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE, TRANS */
/* INTEGER INFO, K, LDA, LDC, LWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), C( LDC, * ), TAU( * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CUNMQR overwrites the general complex M-by-N matrix C with */
/* > */
/* > SIDE = 'L' SIDE = 'R' */
/* > TRANS = 'N': Q * C C * Q */
/* > TRANS = 'C': Q**H * C C * Q**H */
/* > */
/* > where Q is a complex unitary matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by CGEQRF. Q is of order M if SIDE = 'L' and of order N */
/* > if SIDE = 'R'. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': apply Q or Q**H from the Left;
 */
/* > = 'R': apply Q or Q**H from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N': No transpose, apply Q;
 */
/* > = 'C': Conjugate transpose, apply Q**H. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of elementary reflectors whose product defines */
/* > the matrix Q. */
/* > If SIDE = 'L', M >= K >= 0;
 */
/* > if SIDE = 'R', N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,K) */
/* > The i-th column must contain the vector which defines the */
/* > elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* > CGEQRF in the first k columns of its array argument A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. */
/* > If SIDE = 'L', LDA >= fla_max(1,M);
 */
/* > if SIDE = 'R', LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is COMPLEX array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by CGEQRF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX array, dimension (LDC,N) */
/* > On entry, the M-by-N matrix C. */
/* > On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If SIDE = 'L', LWORK >= fla_max(1,N);
 */
/* > if SIDE = 'R', LWORK >= fla_max(1,M). */
/* > For good performance, LWORK should generally be larger. */
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
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
void cunmqr_(char *side, char *trans, integer *m, integer *n, integer *k, complex *a, integer *lda,
             complex *tau, complex *c__, integer *ldc, complex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256,
             "cunmqr inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS
             ", lda %" FLA_IS ", ldc %" FLA_IS ", lwork %" FLA_IS "",
             *side, *trans, *m, *n, *k, *lda, *ldc, *lwork);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    extern void fla_cunmqr(char *side, char *trans, integer *m, integer *n, integer *k, complex *a,
                           integer *lda, complex *tau, complex *c__, integer *ldc, complex *work,
                           integer *lwork, integer *info);

    fla_cunmqr(side, trans, m, n, k, a, lda, tau, c__, ldc, work, lwork, info);
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
}

void fla_cunmqr(char *side, char *trans, integer *m, integer *n, integer *k, complex *a,
                integer *lda, complex *tau, complex *c__, integer *ldc, complex *work,
                integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__4, i__5;
    char ch__1[2];
    /* Builtin functions */
    /* Subroutine */

    /* Local variables */
    integer i__, i1, i2, i3, ib, ic, jc, nb, mi, ni, nq, nw, iwt;
    logical left;
    extern logical lsame_(char *, char *, integer, integer);
    integer nbmin, iinfo;
    extern /* Subroutine */
        void
        cunm2r_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *,
                complex *, integer *, complex *, integer *),
        clarfb_(char *, char *, char *, char *, integer *, integer *, integer *, complex *,
                integer *, complex *, integer *, complex *, integer *, complex *, integer *),
        clarft_(char *, char *, integer *, integer *, complex *, integer *, complex *, complex *,
                integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern real sroundup_lwork(integer *);
    logical notran;
    integer ldwork, lwkopt;
    logical lquery;
#if FLA_OPENMP_MULTITHREADING
    int thread_id, actual_num_threads;
    integer index, mi_sub, ni_sub;
#endif
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */

    /* Initialize global context data */
    aocl_fla_init();

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    *info = 0;
    left = lsame_(side, "L", 1, 1);
    notran = lsame_(trans, "N", 1, 1);
    lquery = *lwork == -1;
    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if(left)
    {
        nq = *m;
        nw = *n;
    }
    else
    {
        nq = *n;
        nw = *m;
    }
    if(!left && !lsame_(side, "R", 1, 1))
    {
        *info = -1;
    }
    else if(!notran && !lsame_(trans, "C", 1, 1))
    {
        *info = -2;
    }
    else if(*m < 0)
    {
        *info = -3;
    }
    else if(*n < 0)
    {
        *info = -4;
    }
    else if(*k < 0 || *k > nq)
    {
        *info = -5;
    }
    else if(*lda < fla_max(1, nq))
    {
        *info = -7;
    }
    else if(*ldc < fla_max(1, *m))
    {
        *info = -10;
    }
    else if(*lwork < fla_max(1, nw) && !lquery)
    {
        *info = -12;
    }
    if(*info == 0)
    {
        /* Compute the workspace requirements */
        /* Computing MIN */
        i__1 = 64;
        i__2 = ilaenv_(&c__1, "CUNMQR", ch__1, m, n, k, &c_n1); // , expr subst
        nb = fla_min(i__1, i__2);
        lwkopt = fla_max(1, nw) * nb + 4160;
        work[1].r = sroundup_lwork(&lwkopt);
        work[1].i = 0.f; // , expr subst
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CUNMQR", &i__1, (ftnlen)6);
        return;
    }
    else if(lquery)
    {
        return;
    }
    /* Quick return if possible */
    if(*m == 0 || *n == 0 || *k == 0)
    {
        work[1].r = 1.f;
        work[1].i = 0.f; // , expr subst
        return;
    }
    nbmin = 2;
    ldwork = nw;
    if(nb > 1 && nb < *k)
    {
        if(*lwork < nw * nb + 4160)
        {
            nb = (*lwork - 4160) / ldwork;
            /* Computing MAX */
            i__1 = 2;
            i__2 = ilaenv_(&c__2, "CUNMQR", ch__1, m, n, k, &c_n1); // , expr subst
            nbmin = fla_max(i__1, i__2);
        }
    }
    if(nb < nbmin || nb >= *k)
    {
        /* Use unblocked code */
        cunm2r_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[c_offset], ldc, &work[1],
                &iinfo);
    }
    else
    {
        /* Use blocked code */
        iwt = nw * nb + 1;
        if(left && !notran || !left && notran)
        {
            i1 = 1;
            i2 = *k;
            i3 = nb;
        }
        else
        {
            i1 = (*k - 1) / nb * nb + 1;
            i2 = 1;
            i3 = -nb;
        }
        if(left)
        {
            ni = *n;
            jc = 1;
        }
        else
        {
            mi = *m;
            ic = 1;
        }
        i__1 = i2;
        i__2 = i3;

#ifdef FLA_OPENMP_MULTITHREADING
        /* Get optimum thread number for DORMLQ*/
        FLA_Thread_optimum(FLA_ORMQR, &actual_num_threads);
#pragma omp parallel num_threads(actual_num_threads) private(i__, thread_id, mi_sub, ni_sub, index)
        {
            thread_id = omp_get_thread_num();
#else
        {
#endif
            for(i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
            {
                /* Computing MIN */
#ifdef FLA_OPENMP_MULTITHREADING
/* Compute triangular factor of the block reflector in a single thread */
#pragma omp single
#endif
                {
                    i__4 = nb;
                    i__5 = *k - i__ + 1; // , expr subst
                    ib = fla_min(i__4, i__5);
                    /* Form the triangular factor of the block reflector */
                    /* H = H(i) H(i+1) . . . H(i+ib-1) */
                    i__4 = nq - i__ + 1;
                    clarft_("Forward", "Columnwise", &i__4, &ib, &a[i__ + i__ * a_dim1], lda,
                            &tau[i__], &work[iwt], &c__65);
                }

                if(left)
                {
                    /* H or H**H is applied to C(i:m,1:n) */
                    mi = *m - i__ + 1;
                    ic = i__;
#ifdef FLA_OPENMP_MULTITHREADING
                    /* Determine the sub partition range of current thread */
                    FLA_Thread_get_subrange(thread_id, actual_num_threads, ni, &ni_sub, &index);
                    mi_sub = mi;
#endif
                }
                else
                {
                    /* H or H**H is applied to C(1:m,i:n) */
                    ni = *n - i__ + 1;
                    jc = i__;
#ifdef FLA_OPENMP_MULTITHREADING
                    /* Determine the sub partition range of current thread */
                    FLA_Thread_get_subrange(thread_id, actual_num_threads, mi, &mi_sub, &index);
                    ni_sub = ni;
#endif
                }
                /* Apply H or H**H */
#ifdef FLA_OPENMP_MULTITHREADING
                if(left)
                    clarfb_(side, trans, "Forward", "Columnwise", &mi_sub, &ni_sub, &ib,
                            &a[i__ + i__ * a_dim1], lda, &work[iwt], &c__65,
                            &c__[ic + (index + jc) * c_dim1], ldc, &work[1 + index], &ldwork);
                else
                    clarfb_(side, trans, "Forward", "Columnwise", &mi_sub, &ni_sub, &ib,
                            &a[i__ + i__ * a_dim1], lda, &work[iwt], &c__65,
                            &c__[index + ic + jc * c_dim1], ldc, &work[1 + index], &ldwork);
#pragma omp barrier
#else
                clarfb_(side, trans, "Forward", "Columnwise", &mi, &ni, &ib, &a[i__ + i__ * a_dim1],
                        lda, &work[iwt], &c__65, &c__[ic + jc * c_dim1], ldc, &work[1], &ldwork);
#endif
                /* L10: */
            }
        }
    }
    work[1].r = sroundup_lwork(&lwkopt);
    work[1].i = 0.f; // , expr subst
    return;
    /* End of CUNMQR */
}
/* cunmqr_ */
