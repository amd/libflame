/* ./cgelqt3.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 = {1.f, 0.f};
/* > \brief \b CGELQT3 */
/* Definition: */
/* =========== */
/* SUBROUTINE CGELQT3( M, N, A, LDA, T, LDT, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, M, N, LDT */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), T( LDT, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGELQT3 recursively computes a LQ factorization of a complex M-by-N */
/* > matrix A, using the compact WY representation of Q. */
/* > */
/* > Based on the algorithm of Elmroth and Gustavson, */
/* > IBM J. Res. Develop. Vol 44 No. 4 July 2000. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M =< N. */
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
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the complex M-by-N matrix A. On exit, the elements on and */
/* > below the diagonal contain the N-by-N lower triangular matrix L;
the */
/* > elements above the diagonal are the rows of V. See below for */
/* > further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is COMPLEX array, dimension (LDT,N) */
/* > The N-by-N upper triangular factor of the block reflector. */
/* > The elements on and above the diagonal contain the block */
/* > reflector T;
the elements below the diagonal are not used. */
/* > See below for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= fla_max(1,N). */
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
/* > \ingroup gelqt3 */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix V stores the elementary reflectors H(i) in the i-th row */
/* > above the diagonal. For example, if M=5 and N=3, the matrix V is */
/* > */
/* > V = ( 1 v1 v1 v1 v1 ) */
/* > ( 1 v2 v2 v2 ) */
/* > ( 1 v3 v3 v3 ) */
/* > */
/* > */
/* > where the vi's represent the vectors which define H(i), which are returned */
/* > in the matrix A. The 1's along the diagonal of V are not stored in A. The */
/* > block reflector H is then given by */
/* > */
/* > H = I - V * T * V**T */
/* > */
/* > where V**T is the transpose of V. */
/* > */
/* > For details of the algorithm, see Elmroth and Gustavson (cited above). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void cgelqt3_(integer *m, integer *n, complex *a, integer *lda, complex *t, integer *ldt,
              integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "cgelqt3 inputs: m %lld, n %lld, lda %lld, ldt %lld", *m, *n, *lda, *ldt);
#else
    snprintf(buffer, 256, "cgelqt3 inputs: m %d, n %d, lda %d, ldt %d", *m, *n, *lda, *ldt);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1;
    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    /* Local variables */
    integer i__, j, i1, j1, m1, m2;
    extern /* Subroutine */
        void
        cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *,
               complex *, integer *, complex *, complex *, integer *);
    integer iinfo;
    extern /* Subroutine */
        void
        ctrmm_(char *, char *, char *, char *, integer *, integer *, complex *, complex *,
               integer *, complex *, integer *),
        clarfg_(integer *, complex *, complex *, integer *, complex *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    /* Function Body */
    *info = 0;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < *m)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -4;
    }
    else if(*ldt < fla_max(1, *m))
    {
        *info = -6;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGELQT3", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(*m == 1)
    {
        /* Compute Householder transform when M=1 */
        clarfg_(n, &a[a_offset], &a[fla_min(2, *n) * a_dim1 + 1], lda, &t[t_offset]);
        i__1 = t_dim1 + 1;
        r_cnjg(&q__1, &t[t_dim1 + 1]);
        t[i__1].r = q__1.r;
        t[i__1].i = q__1.i; // , expr subst
    }
    else
    {
        /* Otherwise, split A into blocks... */
        m1 = *m / 2;
        m2 = *m - m1;
        /* Computing MIN */
        i__1 = m1 + 1;
        i1 = fla_min(i__1, *m);
        /* Computing MIN */
        i__1 = *m + 1;
        j1 = fla_min(i__1, *n);
        /* Compute A(1:M1,1:N) <- (Y1,R1,T1), where Q1 = I - Y1 T1 Y1^H */
        cgelqt3_(&m1, n, &a[a_offset], lda, &t[t_offset], ldt, &iinfo);
        /* Compute A(J1:M,1:N) = A(J1:M,1:N) Q1^H [workspace: T(1:N1,J1:N)] */
        i__1 = m2;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i__2 = m1;
            for(j = 1; j <= i__2; ++j)
            {
                i__3 = i__ + m1 + j * t_dim1;
                i__4 = i__ + m1 + j * a_dim1;
                t[i__3].r = a[i__4].r;
                t[i__3].i = a[i__4].i; // , expr subst
            }
        }
        ctrmm_("R", "U", "C", "U", &m2, &m1, &c_b1, &a[a_offset], lda, &t[i1 + t_dim1], ldt);
        i__1 = *n - m1;
        cgemm_("N", "C", &m2, &m1, &i__1, &c_b1, &a[i1 + i1 * a_dim1], lda, &a[i1 * a_dim1 + 1],
               lda, &c_b1, &t[i1 + t_dim1], ldt);
        ctrmm_("R", "U", "N", "N", &m2, &m1, &c_b1, &t[t_offset], ldt, &t[i1 + t_dim1], ldt);
        i__1 = *n - m1;
        q__1.r = -1.f;
        q__1.i = -0.f; // , expr subst
        cgemm_("N", "N", &m2, &i__1, &m1, &q__1, &t[i1 + t_dim1], ldt, &a[i1 * a_dim1 + 1], lda,
               &c_b1, &a[i1 + i1 * a_dim1], lda);
        ctrmm_("R", "U", "N", "U", &m2, &m1, &c_b1, &a[a_offset], lda, &t[i1 + t_dim1], ldt);
        i__1 = m2;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i__2 = m1;
            for(j = 1; j <= i__2; ++j)
            {
                i__3 = i__ + m1 + j * a_dim1;
                i__4 = i__ + m1 + j * a_dim1;
                i__5 = i__ + m1 + j * t_dim1;
                q__1.r = a[i__4].r - t[i__5].r;
                q__1.i = a[i__4].i - t[i__5].i; // , expr subst
                a[i__3].r = q__1.r;
                a[i__3].i = q__1.i; // , expr subst
                i__3 = i__ + m1 + j * t_dim1;
                t[i__3].r = 0.f;
                t[i__3].i = 0.f; // , expr subst
            }
        }
        /* Compute A(J1:M,J1:N) <- (Y2,R2,T2) where Q2 = I - Y2 T2 Y2^H */
        i__1 = *n - m1;
        cgelqt3_(&m2, &i__1, &a[i1 + i1 * a_dim1], lda, &t[i1 + i1 * t_dim1], ldt, &iinfo);
        /* Compute T3 = T(J1:N1,1:N) = -T1 Y1^H Y2 T2 */
        i__1 = m2;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i__2 = m1;
            for(j = 1; j <= i__2; ++j)
            {
                i__3 = j + (i__ + m1) * t_dim1;
                i__4 = j + (i__ + m1) * a_dim1;
                t[i__3].r = a[i__4].r;
                t[i__3].i = a[i__4].i; // , expr subst
            }
        }
        ctrmm_("R", "U", "C", "U", &m1, &m2, &c_b1, &a[i1 + i1 * a_dim1], lda, &t[i1 * t_dim1 + 1],
               ldt);
        i__1 = *n - *m;
        cgemm_("N", "C", &m1, &m2, &i__1, &c_b1, &a[j1 * a_dim1 + 1], lda, &a[i1 + j1 * a_dim1],
               lda, &c_b1, &t[i1 * t_dim1 + 1], ldt);
        q__1.r = -1.f;
        q__1.i = -0.f; // , expr subst
        ctrmm_("L", "U", "N", "N", &m1, &m2, &q__1, &t[t_offset], ldt, &t[i1 * t_dim1 + 1], ldt);
        ctrmm_("R", "U", "N", "N", &m1, &m2, &c_b1, &t[i1 + i1 * t_dim1], ldt, &t[i1 * t_dim1 + 1],
               ldt);
        /* Y = (Y1,Y2);
        L = [ L1 0 ];
        T = [T1 T3] */
        /* [ A(1:N1,J1:N) L2 ] [ 0 T2] */
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGELQT3 */
}
/* cgelqt3_ */
