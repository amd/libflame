/* ../netlib/clag2z.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
/******************************************************************************
 * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#include "FLA_f2c.h"
/* > \brief \b CLAG2Z converts a scomplex single precision matrix to a scomplex double precision
 * matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAG2Z + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clag2z.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clag2z.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clag2z.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAG2Z( M, N, SA, LDSA, A, LDA, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDSA, M, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX SA( LDSA, * ) */
/* COMPLEX*16 A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAG2Z converts a COMPLEX matrix, SA, to a COMPLEX*16 matrix, A. */
/* > */
/* > Note that while it is possible to overflow while converting */
/* > from double to single, it is not possible to overflow when */
/* > converting from single to double. */
/* > */
/* > This is an auxiliary routine so there is no argument checking. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of lines of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] SA */
/* > \verbatim */
/* > SA is COMPLEX array, dimension (LDSA,N) */
/* > On entry, the M-by-N coefficient matrix SA. */
/* > \endverbatim */
/* > */
/* > \param[in] LDSA */
/* > \verbatim */
/* > LDSA is INTEGER */
/* > The leading dimension of the array SA. LDSA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On exit, the M-by-N coefficient matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void clag2z_(aocl_int_t *m, aocl_int_t *n, scomplex *sa, aocl_int_t *ldsa, dcomplex *a,
             aocl_int_t *lda, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clag2z(m, n, sa, ldsa, a, lda, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldsa_64 = *ldsa;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t info_64 = *info;

    aocl_lapack_clag2z(&m_64, &n_64, sa, &ldsa_64, a, &lda_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_clag2z(aocl_int64_t *m, aocl_int64_t *n, scomplex *sa, aocl_int64_t *ldsa,
                        dcomplex *a, aocl_int64_t *lda, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clag2z inputs: m %lld, n %lld, ldsa %lld, lda %lld", *m, *n, *ldsa,
             *lda);
#else
    snprintf(buffer, 256, "clag2z inputs: m %d, n %d, ldsa %d, lda %d", *m, *n, *ldsa, *lda);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t sa_dim1, sa_offset, a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    aocl_int64_t i__, j;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    sa_dim1 = *ldsa;
    sa_offset = 1 + sa_dim1;
    sa -= sa_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    /* Function Body */
    *info = 0;
    i__1 = *n;
    for(j = 1; j <= i__1; ++j)
    {
        i__2 = *m;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            i__3 = i__ + j * a_dim1;
            i__4 = i__ + j * sa_dim1;
            a[i__3].real = sa[i__4].real;
            a[i__3].imag = sa[i__4].imag; // , expr subst
            /* L10: */
        }
        /* L20: */
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLAG2Z */
}
/* clag2z_ */
