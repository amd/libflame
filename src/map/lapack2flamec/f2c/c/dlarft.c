/* ../netlib/dlarft.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
/******************************************************************************
 * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#include "FLA_f2c.h" /* Table of constant values */
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif

static integer c__1 = 1;
static doublereal c_b6 = 1.;
static doublereal c_b0 = 0.;
/* > \brief \b DLARFT forms the triangular factor T of a block reflector H = I - vtvH */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLARFT + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarft.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarft.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarft.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIRECT, STOREV */
/* INTEGER K, LDT, LDV, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION T( LDT, * ), TAU( * ), V( LDV, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARFT forms the triangular factor T of a real block reflector H */
/* > of order n, which is defined as a product of k elementary reflectors. */
/* > */
/* > If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
 */
/* > */
/* > If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular. */
/* > */
/* > If STOREV = 'C', the vector which defines the elementary reflector */
/* > H(i) is stored in the i-th column of the array V, and */
/* > */
/* > H = I - V * T * V**T */
/* > */
/* > If STOREV = 'R', the vector which defines the elementary reflector */
/* > H(i) is stored in the i-th row of the array V, and */
/* > */
/* > H = I - V**T * T * V */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] DIRECT */
/* > \verbatim */
/* > DIRECT is CHARACTER*1 */
/* > Specifies the order in which the elementary reflectors are */
/* > multiplied to form the block reflector: */
/* > = 'F': H = H(1) H(2) . . . H(k) (Forward) */
/* > = 'B': H = H(k) . . . H(2) H(1) (Backward) */
/* > \endverbatim */
/* > */
/* > \param[in] STOREV */
/* > \verbatim */
/* > STOREV is CHARACTER*1 */
/* > Specifies how the vectors which define the elementary */
/* > reflectors are stored (see also Further Details): */
/* > = 'C': columnwise */
/* > = 'R': rowwise */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the block reflector H. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The order of the triangular factor T (= the number of */
/* > elementary reflectors). K >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* > V is DOUBLE PRECISION array, dimension */
/* > (LDV,K) if STOREV = 'C' */
/* > (LDV,N) if STOREV = 'R' */
/* > The matrix V. See further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > The leading dimension of the array V. */
/* > If STOREV = 'C', LDV >= fla_max(1,N);
if STOREV = 'R', LDV >= K. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is DOUBLE PRECISION array, dimension (LDT,K) */
/* > The k by k triangular factor T of the block reflector. */
/* > If DIRECT = 'F', T is upper triangular;
if DIRECT = 'B', T is */
/* > lower triangular. The rest of the array is not used. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= K. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup doubleOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The shape of the matrix V and the storage of the vectors which define */
/* > the H(i) is best illustrated by the following example with n = 5 and */
/* > k = 3. The elements equal to 1 are not stored. */
/* > */
/* > DIRECT = 'F' and STOREV = 'C': DIRECT = 'F' and STOREV = 'R': */
/* > */
/* > V = ( 1 ) V = ( 1 v1 v1 v1 v1 ) */
/* > ( v1 1 ) ( 1 v2 v2 v2 ) */
/* > ( v1 v2 1 ) ( 1 v3 v3 ) */
/* > ( v1 v2 v3 ) */
/* > ( v1 v2 v3 ) */
/* > */
/* > DIRECT = 'B' and STOREV = 'C': DIRECT = 'B' and STOREV = 'R': */
/* > */
/* > V = ( v1 v2 v3 ) V = ( v1 v1 1 ) */
/* > ( v1 v2 v3 ) ( v2 v2 v2 1 ) */
/* > ( 1 v2 v3 ) ( v3 v3 v3 v3 1 ) */
/* > ( 1 v3 ) */
/* > ( 1 ) */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void dlarft_(char *direct, char *storev, integer *n, integer *k, doublereal *v, integer *ldv,
             doublereal *tau, doublereal *t, integer *ldt)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlarft inputs: direct %c, storev %c, n %" FLA_IS ", k %" FLA_IS
                      ", ldv %" FLA_IS ", ldt %" FLA_IS "",
                      *direct, *storev, *n, *k, *ldv, *ldt);
    /* System generated locals */
    integer t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3;
    doublereal d__1;
    /* Local variables */
    integer i__, j, prevlastv;
    integer lastv;
#if !FLA_ENABLE_AOCL_BLAS
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *,
               integer *, doublereal *, doublereal *, integer *);
    extern /* Subroutine */
        void
        dtrmv_(char *, char *, char *, integer *, doublereal *, integer *, doublereal *, integer *);
#ifdef FLA_ENABLE_AMD_OPT
    extern void dtrmm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                       doublereal *, integer *, doublereal *, integer *);
    extern void dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                       integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
#endif
#endif

    /* -- LAPACK auxiliary routine (version 3.7.0) -- */
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
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --tau;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    /* Function Body */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#ifdef FLA_ENABLE_AMD_OPT
	aocl_fla_init();
#endif
    if(lsame_(direct, "F", 1, 1))
    {
        prevlastv = *n;
        i__1 = *k;
#ifdef FLA_ENABLE_AMD_OPT
        if(lsame_(storev, "C", 1, 1))
        {
            integer nb = FLA_DLARFT_BLOCK_NB;
            integer f_ncols = fla_min(nb, *k);
            doublereal n_tau;
            doublereal diag_elem;
            integer block_last_i;

            /*
             *    Let V be the matrix (n x k) of the elementary reflectors
             *
             *
             *    V is partitioned as follows:
             *
             *    V = || V11  V12 ||
             *        || V21  V22 ||
             *        || V31  V32 ||
             *
             *    Where,
             *         V11 is unit lower triangular
             *         V21 is rectangular
             *         V31 is rectangular
             *         V12 is rectangular
             *         V22 is rectangular
             *         V32 is rectangular
             *
             *    Let T be the matrix (k x k) of the elementary reflectors
             *
             *    T is paritioned as follows:
             *
             *    T = || T11  T12 ||
             *        || T21  T22 ||
             *
             *    Where,
             *         T11 is upper triangular
             *         T21 is zero
             *         T12 is rectangular
             *         T22 is upper triangular
             *
             *
             */

            /* For the first column of T, we have T11(1,1) = tau[1] */
            t[t_dim1 + 1] = tau[1];

            /*
             * Generate T11 using unblocked algorithm
             *
             */
            for(i__ = 2; i__ <= f_ncols; ++i__)
            {
                prevlastv = fla_max(prevlastv, i__);
                integer m_lower_triangular = prevlastv - i__ + 1;
                integer n_left = i__ - 1;
                n_tau = -tau[i__];
                diag_elem = v[i__ + i__ * v_dim1];
                /* Explicitly make the diagonal element of V11 equal to 1 */
                v[i__ + i__ * v_dim1] = 1.;
                /* T11(1:i-1,i) := - tau(i) * V11(i:j,1:i-1)**T * V11(i:j,i) */
#if FLA_ENABLE_AOCL_BLAS
                if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
                {
                    bli_dgemv_t_zen4_int(BLIS_CONJUGATE, BLIS_NO_CONJUGATE, m_lower_triangular,
                                         n_left, &n_tau, &v[i__ + v_dim1], 1, *ldv,
                                         &v[i__ + i__ * v_dim1], c__1,
                                         &c_b0, &t[i__ * t_dim1 + 1], c__1, NULL);
                }
                else
#endif
                {
                    dgemv_("Transpose", &m_lower_triangular, &n_left, &n_tau, &v[i__ + v_dim1], ldv,
                           &v[i__ + i__ * v_dim1], &c__1, &c_b0, &t[i__ * t_dim1 + 1], &c__1);
                }
                /* Restore the diagonal element of V11 */
                v[i__ + i__ * v_dim1] = diag_elem;
                /* T11(1:i-1,i) := T11(1:i-1,1:i-1) * T11(1:i-1,i) */
                dtrmv_("Upper", "No transpose", "Non-unit", &n_left, &t[t_offset], ldt,
                       &t[i__ * t_dim1 + 1], &c__1);
                /* T11(i,i) := tau(i) */
                t[i__ + i__ * t_dim1] = tau[i__];
            }

            /* Process the remaining blocks from column nb + 1 to k */
            for(i__ = nb + 1; i__ <= *k; i__ += nb)
            {
                /* Using gemm for partial update of T12 */
                block_last_i = fla_min(i__ + nb - 1, *k);
                integer n_v32 = block_last_i - i__ + 1;
                integer m_v31 = fla_max(*n, block_last_i) - block_last_i;
                integer n_v31 = i__ - 1;

                /* T12 = V31**T * V32 */
                dgemm_("Transpose", "No transpose", &n_v31, &n_v32, &m_v31, &c_b6,
                       &v[block_last_i + 1 + v_dim1], ldv, &v[block_last_i + 1 + i__ * v_dim1], ldv,
                       &c_b0, &t[i__ * t_dim1 + 1], ldt);

                for(j = i__; j <= block_last_i; ++j)
                {
                    integer m_v22_j = block_last_i - j + 1;
                    n_tau = -tau[j];
                    diag_elem = v[j + j * v_dim1];
                    /* Explicitly make the diagonal element of V22 equal to 1 */
                    v[j + j * v_dim1] = 1.;
                    /* Update T12
                       T12(:, j) = -tau[j] * T12(:, j) -tau[j] * V21**T * V22(:, j)
                     */
#if FLA_ENABLE_AOCL_BLAS
                    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
                    {
                        double n_tau_d = n_tau;
                        bli_dgemv_t_zen4_int(BLIS_CONJUGATE, BLIS_NO_CONJUGATE, m_v22_j, n_v31, &n_tau,
                                             &v[j + v_dim1], 1, *ldv,
                                             &v[j + j * v_dim1], c__1,
                                             &n_tau_d, &t[j * t_dim1 + 1], c__1, NULL);
                    }
                    else
#endif
                    {
                        dgemv_("Transpose", &m_v22_j, &n_v31, &n_tau, &v[j + v_dim1], ldv,
                               &v[j + j * v_dim1], &c__1, &n_tau, &t[j * t_dim1 + 1], &c__1);
                    }

                    /* V22_32 = || V22 ||
                     *          || V32 ||
                     */

                    integer m_v22_32_j = fla_max(*n, j) - j + 1;
                    integer n_v22_32_j = j - i__;

                    /* Update T22
                     * T22(:, j) =  -tau[j] * V22_32(:,1:j-1)**T * V22_32(:, j)
                     */
#if FLA_ENABLE_AOCL_BLAS
                    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
                    {
                        bli_dgemv_t_zen4_int(BLIS_CONJUGATE, BLIS_NO_CONJUGATE, m_v22_32_j,
                                             n_v22_32_j, &n_tau, &v[j + v_dim1 * i__], 1, *ldv,
                                             &v[j + j * v_dim1], c__1,
                                             &c_b0, &t[j * t_dim1 + i__], c__1, NULL);
                    }
                    else
#endif
                    {
                        dgemv_("Transpose", &m_v22_32_j, &n_v22_32_j, &n_tau, &v[j + v_dim1 * i__], ldv,
                               &v[j + j * v_dim1], &c__1, &c_b0, &t[j * t_dim1 + i__], &c__1);
                    }

                    /* Restore the diagonal element of V22 */
                    v[j + j * v_dim1] = diag_elem;
                    /* T22(j, j) = tau[j] */
                    t[j + j * t_dim1] = tau[j];
                }

                integer m_t12 = i__ - 1;
                integer n_t12 = block_last_i - i__ + 1;

                /* Update T12
                 * T12 = T11 * T12
                 */
                dtrmm_("Left", "Upper", "No transpose", "Non-unit", &m_t12, &n_t12, &c_b6,
                       &t[t_offset], ldt, &t[i__ * t_dim1 + 1], ldt);

                for(j = i__ + 1; j <= block_last_i; ++j)
                {
                    integer n_t12_j = j - i__;
                    /* Update T12
                     * T12(:, j) = T12(:,j) + T12 * T12(:, j)
                     */
#if FLA_ENABLE_AOCL_BLAS
                    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
                    {
                        bli_dgemv_n_zen4_int_40x2_st(BLIS_NO_TRANSPOSE, BLIS_NO_CONJUGATE, m_t12,
                                                     n_t12_j, &c_b6, &t[1 + i__ * t_dim1], 1, *ldt,
                                                     &t[i__ + j * t_dim1], c__1,
                                                     &c_b6, &t[j * t_dim1 + 1], c__1, NULL);
                    }
                    else
#endif
                    {
                        dgemv_("No transpose", &m_t12, &n_t12_j, &c_b6, &t[1 + i__ * t_dim1], ldt,
                               &t[i__ + j * t_dim1], &c__1, &c_b6, &t[j * t_dim1 + 1], &c__1);
                    }
                    /*
                     * Update T22
                     * T22(:, j) = T22(1:j-1, j) * T22(:, j)
                     */
                    dtrmv_("Upper", "No transpose", "Non-unit", &n_t12_j, &t[i__ + i__ * t_dim1],
                           ldt, &t[i__ + j * t_dim1], &c__1);
                }
            }
        }
        else
        {
#endif
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                prevlastv = fla_max(i__, prevlastv);
                if(tau[i__] == 0.)
                {
                    /* H(i) = I */
                    i__2 = i__;
                    for(j = 1; j <= i__2; ++j)
                    {
                        t[j + i__ * t_dim1] = 0.;
                    }
                }
                else
                {
                    /* general case */
                    if(lsame_(storev, "C", 1, 1))
                    {
                        /* Skip any trailing zeros. */
                        i__2 = i__ + 1;
                        for(lastv = *n; lastv >= i__2; --lastv)
                        {
                            if(v[lastv + i__ * v_dim1] != 0.)
                            {
                                break;
                            }
                        }
                        i__2 = i__ - 1;
                        for(j = 1; j <= i__2; ++j)
                        {
                            t[j + i__ * t_dim1] = -tau[i__] * v[i__ + j * v_dim1];
                        }
                        j = fla_min(lastv, prevlastv);
                        /* T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i) */
                        i__2 = j - i__;
                        i__3 = i__ - 1;
                        d__1 = -tau[i__];
#if FLA_ENABLE_AOCL_BLAS
                        if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
                        {
                            bli_dgemv_t_zen4_int(BLIS_CONJUGATE, BLIS_NO_CONJUGATE, i__2, i__3, &d__1,
                                                 &v[i__ + 1 + v_dim1], 1, *ldv,
                                                 &v[i__ + 1 + i__ * v_dim1], c__1,
                                                 &c_b6, &t[i__ * t_dim1 + 1], c__1, NULL);
                        }
                        else
#endif
                        {
                            dgemv_("Transpose", &i__2, &i__3, &d__1, &v[i__ + 1 + v_dim1], ldv,
                                   &v[i__ + 1 + i__ * v_dim1], &c__1, &c_b6, &t[i__ * t_dim1 + 1],
                                &c__1);
                        }
                    }
                    else
                    {
                        /* Skip any trailing zeros. */
                        i__2 = i__ + 1;
                        for(lastv = *n; lastv >= i__2; --lastv)
                        {
                            if(v[i__ + lastv * v_dim1] != 0.)
                            {
                                break;
                            }
                        }
                        i__2 = i__ - 1;
                        for(j = 1; j <= i__2; ++j)
                        {
                            t[j + i__ * t_dim1] = -tau[i__] * v[j + i__ * v_dim1];
                        }
                        j = fla_min(lastv, prevlastv);
                        /* T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**T */
                        i__2 = i__ - 1;
                        i__3 = j - i__;
                        d__1 = -tau[i__];
#if FLA_ENABLE_AOCL_BLAS
                        if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512) && *ldv > 0)
                        {
                            bli_dgemv_n_zen4_int_40x2_st(BLIS_NO_TRANSPOSE, BLIS_NO_CONJUGATE, i__2, i__3, &d__1,
                                                 &v[(i__ + 1) * v_dim1 + 1], 1, *ldv,
                                                 &v[i__ + (i__ + 1) * v_dim1], *ldv,
                                                 &c_b6, &t[i__ * t_dim1 + 1], c__1, NULL);
                        }
                        else
#endif
                        {
                            dgemv_("No transpose", &i__2, &i__3, &d__1, &v[(i__ + 1) * v_dim1 + 1], ldv,
                                   &v[i__ + (i__ + 1) * v_dim1], ldv, &c_b6, &t[i__ * t_dim1 + 1],
                                   &c__1);
                        }
                    }
                    /* T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i) */
                    i__2 = i__ - 1;
                    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[t_offset], ldt,
                           &t[i__ * t_dim1 + 1], &c__1);
                    t[i__ + i__ * t_dim1] = tau[i__];
                    if(i__ > 1)
                    {
                        prevlastv = fla_max(prevlastv, lastv);
                    }
                    else
                    {
                        prevlastv = lastv;
                    }
                }
            }
#ifdef FLA_ENABLE_AMD_OPT
        }
#endif
    }
    else
    {
        prevlastv = 1;
        for(i__ = *k; i__ >= 1; --i__)
        {
            if(tau[i__] == 0.)
            {
                /* H(i) = I */
                i__1 = *k;
                for(j = i__; j <= i__1; ++j)
                {
                    t[j + i__ * t_dim1] = 0.;
                }
            }
            else
            {
                /* general case */
                if(i__ < *k)
                {
                    if(lsame_(storev, "C", 1, 1))
                    {
                        /* Skip any leading zeros. */
                        i__1 = i__ - 1;
                        for(lastv = 1; lastv <= i__1; ++lastv)
                        {
                            if(v[lastv + i__ * v_dim1] != 0.)
                            {
                                break;
                            }
                        }
                        i__1 = *k;
                        for(j = i__ + 1; j <= i__1; ++j)
                        {
                            t[j + i__ * t_dim1] = -tau[i__] * v[*n - *k + i__ + j * v_dim1];
                        }
                        j = fla_max(lastv, prevlastv);
                        /* T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**T * V(j:n-k+i,i) */
                        i__1 = *n - *k + i__ - j;
                        i__2 = *k - i__;
                        d__1 = -tau[i__];
#if FLA_ENABLE_AOCL_BLAS
                        if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
                        {
                            bli_dgemv_t_zen4_int(BLIS_CONJUGATE, BLIS_NO_CONJUGATE, i__1, i__2, &d__1,
                                                 &v[j + (i__ + 1) * v_dim1], 1, *ldv,
                                                 &v[j + i__ * v_dim1], c__1,
                                                 &c_b6, &t[i__ + 1 + i__ * t_dim1], c__1, NULL);
                        }
                        else
#endif
                        {
                            dgemv_("Transpose", &i__1, &i__2, &d__1, &v[j + (i__ + 1) * v_dim1], ldv,
                                   &v[j + i__ * v_dim1], &c__1, &c_b6, &t[i__ + 1 + i__ * t_dim1],
                                   &c__1);
                        }
                    }
                    else
                    {
                        /* Skip any leading zeros. */
                        i__1 = i__ - 1;
                        for(lastv = 1; lastv <= i__1; ++lastv)
                        {
                            if(v[i__ + lastv * v_dim1] != 0.)
                            {
                                break;
                            }
                        }
                        i__1 = *k;
                        for(j = i__ + 1; j <= i__1; ++j)
                        {
                            t[j + i__ * t_dim1] = -tau[i__] * v[j + (*n - *k + i__) * v_dim1];
                        }
                        j = fla_max(lastv, prevlastv);
                        /* T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**T */
                        i__1 = *k - i__;
                        i__2 = *n - *k + i__ - j;
                        d__1 = -tau[i__];
#if FLA_ENABLE_AOCL_BLAS
                        if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512) && *ldv > 0)
                        {
                            bli_dgemv_n_zen4_int_40x2_st(BLIS_NO_TRANSPOSE, BLIS_NO_CONJUGATE, i__1, i__2, &d__1,
                                                         &v[i__ + 1 + j * v_dim1], 1, *ldv,
                                                         &v[i__ + j * v_dim1], *ldv,
                                                         &c_b6, &t[i__ + 1 + i__ * t_dim1], c__1, NULL);
                        }
                        else
#endif
                        {
                            dgemv_("No transpose", &i__1, &i__2, &d__1, &v[i__ + 1 + j * v_dim1], ldv,
                                   &v[i__ + j * v_dim1], ldv, &c_b6, &t[i__ + 1 + i__ * t_dim1], &c__1);
                        }
                    }
                    /* T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i) */
                    i__1 = *k - i__;
                    dtrmv_("Lower", "No transpose", "Non-unit", &i__1,
                           &t[i__ + 1 + (i__ + 1) * t_dim1], ldt, &t[i__ + 1 + i__ * t_dim1],
                           &c__1);
                    if(i__ > 1)
                    {
                        prevlastv = fla_min(prevlastv, lastv);
                    }
                    else
                    {
                        prevlastv = lastv;
                    }
                }
                t[i__ + i__ * t_dim1] = tau[i__];
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DLARFT */
}
/* dlarft_ */
