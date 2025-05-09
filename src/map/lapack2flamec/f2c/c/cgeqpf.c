/* ../netlib/cgeqpf.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CGEQPF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGEQPF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqpf.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqpf.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqpf.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGEQPF( M, N, A, LDA, JPVT, TAU, WORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER JPVT( * ) */
/* REAL RWORK( * ) */
/* COMPLEX A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine CGEQP3. */
/* > */
/* > CGEQPF computes a QR factorization with column pivoting of a */
/* > complex M-by-N matrix A: A*P = Q*R. */
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
/* > The number of columns of the matrix A. N >= 0 */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, the upper triangle of the array contains the */
/* > fla_min(M,N)-by-N upper triangular matrix R;
the elements */
/* > below the diagonal, together with the array TAU, */
/* > represent the unitary matrix Q as a product of */
/* > fla_min(m,n) elementary reflectors. */
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
/* > On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted */
/* > to the front of A*P (a leading column);
if JPVT(i) = 0, */
/* > the i-th column of A is a free column. */
/* > On exit, if JPVT(i) = k, then the i-th column of A*P */
/* > was the k-th column of A. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX array, dimension (fla_min(M,N)) */
/* > The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (2*N) */
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
/* > \date November 2011 */
/* > \ingroup complexGEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix Q is represented as a product of elementary reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(n) */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar, and v is a complex vector with */
/* > v(1:i-1) = 0 and v(i) = 1;
v(i+1:m) is stored on exit in A(i+1:m,i). */
/* > */
/* > The matrix P is represented in jpvt as follows: If */
/* > jpvt(j) = i */
/* > then the jth column of P is the ith canonical unit vector. */
/* > */
/* > Partial column norm updating strategy modified by */
/* > Z. Drmac and Z. Bujanovic, Dept. of Mathematics, */
/* > University of Zagreb, Croatia. */
/* > -- April 2011 -- */
/* > For more details see LAPACK Working Note 176. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void cgeqpf_(integer *m, integer *n, complex *a, integer *lda, integer *jpvt, complex *tau,
             complex *work, real *rwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256, "cgeqpf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
             *lda);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1, r__2;
    complex q__1;
    /* Builtin functions */
    double sqrt(doublereal);
    void r_cnjg(complex *, complex *);
    double c_abs(complex *);
    /* Local variables */
    integer i__, j, ma, mn;
    complex aii;
    integer pvt;
    real temp, temp2, tol3z;
    extern /* Subroutine */
        void
        clarf_(char *, integer *, integer *, complex *, integer *, complex *, complex *, integer *,
               complex *),
        cswap_(integer *, complex *, integer *, complex *, integer *);
    integer itemp;
    extern /* Subroutine */
        void
        cgeqr2_(integer *, integer *, complex *, integer *, complex *, complex *, integer *);
    extern real scnrm2_(integer *, complex *, integer *);
    extern /* Subroutine */
        void
        cunm2r_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *,
                complex *, integer *, complex *, integer *),
        clarfg_(integer *, complex *, complex *, integer *, complex *);
    extern real slamch_(char *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer isamax_(integer *, real *, integer *);
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
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
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
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
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGEQPF", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    mn = fla_min(*m, *n);
    tol3z = sqrt(slamch_("Epsilon"));
    /* Move initial columns up front */
    itemp = 1;
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if(jpvt[i__] != 0)
        {
            if(i__ != itemp)
            {
                cswap_(m, &a[i__ * a_dim1 + 1], &c__1, &a[itemp * a_dim1 + 1], &c__1);
                jpvt[i__] = jpvt[itemp];
                jpvt[itemp] = i__;
            }
            else
            {
                jpvt[i__] = i__;
            }
            ++itemp;
        }
        else
        {
            jpvt[i__] = i__;
        }
        /* L10: */
    }
    --itemp;
    /* Compute the QR factorization and update remaining columns */
    if(itemp > 0)
    {
        ma = fla_min(itemp, *m);
        cgeqr2_(m, &ma, &a[a_offset], lda, &tau[1], &work[1], info);
        if(ma < *n)
        {
            i__1 = *n - ma;
            cunm2r_("Left", "Conjugate transpose", m, &i__1, &ma, &a[a_offset], lda, &tau[1],
                    &a[(ma + 1) * a_dim1 + 1], lda, &work[1], info);
        }
    }
    if(itemp < mn)
    {
        /* Initialize partial column norms. The first n elements of */
        /* work store the exact column norms. */
        i__1 = *n;
        for(i__ = itemp + 1; i__ <= i__1; ++i__)
        {
            i__2 = *m - itemp;
            rwork[i__] = scnrm2_(&i__2, &a[itemp + 1 + i__ * a_dim1], &c__1);
            rwork[*n + i__] = rwork[i__];
            /* L20: */
        }
        /* Compute factorization */
        i__1 = mn;
        for(i__ = itemp + 1; i__ <= i__1; ++i__)
        {
            /* Determine ith pivot column and swap if necessary */
            i__2 = *n - i__ + 1;
            pvt = i__ - 1 + isamax_(&i__2, &rwork[i__], &c__1);
            if(pvt != i__)
            {
                cswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &c__1);
                itemp = jpvt[pvt];
                jpvt[pvt] = jpvt[i__];
                jpvt[i__] = itemp;
                rwork[pvt] = rwork[i__];
                rwork[*n + pvt] = rwork[*n + i__];
            }
            /* Generate elementary reflector H(i) */
            i__2 = i__ + i__ * a_dim1;
            aii.r = a[i__2].r;
            aii.i = a[i__2].i; // , expr subst
            i__2 = *m - i__ + 1;
            /* Computing MIN */
            i__3 = i__ + 1;
            clarfg_(&i__2, &aii, &a[fla_min(i__3, *m) + i__ * a_dim1], &c__1, &tau[i__]);
            i__2 = i__ + i__ * a_dim1;
            a[i__2].r = aii.r;
            a[i__2].i = aii.i; // , expr subst
            if(i__ < *n)
            {
                /* Apply H(i) to A(i:m,i+1:n) from the left */
                i__2 = i__ + i__ * a_dim1;
                aii.r = a[i__2].r;
                aii.i = a[i__2].i; // , expr subst
                i__2 = i__ + i__ * a_dim1;
                a[i__2].r = 1.f;
                a[i__2].i = 0.f; // , expr subst
                i__2 = *m - i__ + 1;
                i__3 = *n - i__;
                r_cnjg(&q__1, &tau[i__]);
                clarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &q__1,
                       &a[i__ + (i__ + 1) * a_dim1], lda, &work[1]);
                i__2 = i__ + i__ * a_dim1;
                a[i__2].r = aii.r;
                a[i__2].i = aii.i; // , expr subst
            }
            /* Update partial column norms */
            i__2 = *n;
            for(j = i__ + 1; j <= i__2; ++j)
            {
                if(rwork[j] != 0.f)
                {
                    /* NOTE: The following 4 lines follow from the analysis in */
                    /* Lapack Working Note 176. */
                    temp = c_abs(&a[i__ + j * a_dim1]) / rwork[j];
                    /* Computing MAX */
                    r__1 = 0.f;
                    r__2 = (temp + 1.f) * (1.f - temp); // , expr subst
                    temp = fla_max(r__1, r__2);
                    /* Computing 2nd power */
                    r__1 = rwork[j] / rwork[*n + j];
                    temp2 = temp * (r__1 * r__1);
                    if(temp2 <= tol3z)
                    {
                        if(*m - i__ > 0)
                        {
                            i__3 = *m - i__;
                            rwork[j] = scnrm2_(&i__3, &a[i__ + 1 + j * a_dim1], &c__1);
                            rwork[*n + j] = rwork[j];
                        }
                        else
                        {
                            rwork[j] = 0.f;
                            rwork[*n + j] = 0.f;
                        }
                    }
                    else
                    {
                        rwork[j] *= sqrt(temp);
                    }
                }
                /* L30: */
            }
            /* L40: */
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGEQPF */
}
/* cgeqpf_ */
