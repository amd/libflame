/* ../netlib/clatdf.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {1.f, 0.f};
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
static real c_b24 = 1.f;
/* > \brief \b CLATDF uses the LU factorization of the n-by-n matrix computed by sgetc2 and computes
 * a contrib ution to the reciprocal Dif-estimate. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLATDF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatdf.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatdf.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatdf.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, */
/* JPIV ) */
/* .. Scalar Arguments .. */
/* INTEGER IJOB, LDZ, N */
/* REAL RDSCAL, RDSUM */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ), JPIV( * ) */
/* COMPLEX RHS( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLATDF computes the contribution to the reciprocal Dif-estimate */
/* > by solving for x in Z * x = b, where b is chosen such that the norm */
/* > of x is as large as possible. It is assumed that LU decomposition */
/* > of Z has been computed by CGETC2. On entry RHS = f holds the */
/* > contribution from earlier solved sub-systems, and on return RHS = x. */
/* > */
/* > The factorization of Z returned by CGETC2 has the form */
/* > Z = P * L * U * Q, where P and Q are permutation matrices. L is lower */
/* > triangular with unit diagonal elements and U is upper triangular. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] IJOB */
/* > \verbatim */
/* > IJOB is INTEGER */
/* > IJOB = 2: First compute an approximative null-vector e */
/* > of Z using CGECON, e is normalized and solve for */
/* > Zx = +-e - f with the sign giving the greater value of */
/* > 2-norm(x). About 5 times as expensive as Default. */
/* > IJOB .ne. 2: Local look ahead strategy where */
/* > all entries of the r.h.s. b is choosen as either +1 or */
/* > -1. Default. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix Z. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (LDZ, N) */
/* > On entry, the LU part of the factorization of the n-by-n */
/* > matrix Z computed by CGETC2: Z = P * L * U * Q */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDA >= fla_max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] RHS */
/* > \verbatim */
/* > RHS is REAL array, dimension (N). */
/* > On entry, RHS contains contributions from other subsystems. */
/* > On exit, RHS contains the solution of the subsystem with */
/* > entries according to the value of IJOB (see above). */
/* > \endverbatim */
/* > */
/* > \param[in,out] RDSUM */
/* > \verbatim */
/* > RDSUM is REAL */
/* > On entry, the sum of squares of computed contributions to */
/* > the Dif-estimate under computation by CTGSYL, where the */
/* > scaling factor RDSCAL (see below) has been factored out. */
/* > On exit, the corresponding sum of squares updated with the */
/* > contributions from the current sub-system. */
/* > If TRANS = 'T' RDSUM is not touched. */
/* > NOTE: RDSUM only makes sense when CTGSY2 is called by CTGSYL. */
/* > \endverbatim */
/* > */
/* > \param[in,out] RDSCAL */
/* > \verbatim */
/* > RDSCAL is REAL */
/* > On entry, scaling factor used to prevent overflow in RDSUM. */
/* > On exit, RDSCAL is updated w.r.t. the current contributions */
/* > in RDSUM. */
/* > If TRANS = 'T', RDSCAL is not touched. */
/* > NOTE: RDSCAL only makes sense when CTGSY2 is called by */
/* > CTGSYL. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N). */
/* > The pivot indices;
for 1 <= i <= N, row i of the */
/* > matrix has been interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] JPIV */
/* > \verbatim */
/* > JPIV is INTEGER array, dimension (N). */
/* > The pivot indices;
for 1 <= j <= N, column j of the */
/* > matrix has been interchanged with column JPIV(j). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > This routine is a further developed implementation of algorithm */
/* > BSOLVE in [1] using complete pivoting in the LU factorization. */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* > Umea University, S-901 87 Umea, Sweden. */
/* > \par References: */
/* ================ */
/* > */
/* > [1] Bo Kagstrom and Lars Westin, */
/* > Generalized Schur Methods with Condition Estimators for */
/* > Solving the Generalized Sylvester Equation, IEEE Transactions */
/* > on Automatic Control, Vol. 34, No. 7, July 1989, pp 745-751. */
/* > */
/* > [2] Peter Poromaa, */
/* > On Efficient and Robust Estimators for the Separation */
/* > between two Regular Matrix Pairs with Applications in */
/* > Condition Estimation. Report UMINF-95.05, Department of */
/* > Computing Science, Umea University, S-901 87 Umea, Sweden, */
/* > 1995. */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void clatdf_(aocl_int_t *ijob, aocl_int_t *n, scomplex *z__, aocl_int_t *ldz, scomplex *rhs,
             real *rdsum, real *rdscal, aocl_int_t *ipiv, aocl_int_t *jpiv)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clatdf(ijob, n, z__, ldz, rhs, rdsum, rdscal, ipiv, jpiv);
#else
    aocl_int64_t ijob_64 = *ijob;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldz_64 = *ldz;

    aocl_lapack_clatdf(&ijob_64, &n_64, z__, &ldz_64, rhs, rdsum, rdscal, ipiv, jpiv);
#endif
}

void aocl_lapack_clatdf(aocl_int64_t *ijob, aocl_int64_t *n, scomplex *z__, aocl_int64_t *ldz,
                        scomplex *rhs, real *rdsum, real *rdscal, aocl_int_t *ipiv, aocl_int_t *jpiv)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clatdf inputs: ijob %lld, n %lld, ldz %lld", *ijob, *n, *ldz);
#else
    snprintf(buffer, 256, "clatdf inputs: ijob %d, n %d, ldz %d", *ijob, *n, *ldz);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    scomplex q__1, q__2, q__3;
    /* Builtin functions */
    void c_div(scomplex *, scomplex *, scomplex *);
    double c_abs(scomplex *);
    void c_sqrt(scomplex *, scomplex *);
    /* Local variables */
    aocl_int64_t i__, j, k;
    scomplex bm, bp, xm[2], xp[2];
    aocl_int64_t info;
    scomplex temp, work[8];
    real scale;
    scomplex pmone;
    real rtemp, sminu, rwork[2], splus;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --rhs;
    --ipiv;
    --jpiv;
    /* Function Body */
    if(*ijob != 2)
    {
        /* Apply permutations IPIV to RHS */
        i__1 = *n - 1;
        aocl_lapack_claswp(&c__1, &rhs[1], ldz, &c__1, &i__1, &ipiv[1], &c__1);
        /* Solve for L-part choosing RHS either to +1 or -1. */
        q__1.real = -1.f;
        q__1.imag = -0.f; // , expr subst
        pmone.real = q__1.real;
        pmone.imag = q__1.imag; // , expr subst
        i__1 = *n - 1;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j;
            q__1.real = rhs[i__2].real + 1.f;
            q__1.imag = rhs[i__2].imag + 0.f; // , expr subst
            bp.real = q__1.real;
            bp.imag = q__1.imag; // , expr subst
            i__2 = j;
            q__1.real = rhs[i__2].real - 1.f;
            q__1.imag = rhs[i__2].imag - 0.f; // , expr subst
            bm.real = q__1.real;
            bm.imag = q__1.imag; // , expr subst
            splus = 1.f;
            /* Lockahead for L- part RHS(1:N-1) = +-1 */
            /* SPLUS and SMIN computed more efficiently than in BSOLVE[1]. */
            i__2 = *n - j;
            aocl_lapack_cdotc_f2c(&q__1, &i__2, &z__[j + 1 + j * z_dim1], &c__1, &z__[j + 1 + j * z_dim1],
                       &c__1);
            splus += q__1.real;
            i__2 = *n - j;
            aocl_lapack_cdotc_f2c(&q__1, &i__2, &z__[j + 1 + j * z_dim1], &c__1, &rhs[j + 1], &c__1);
            sminu = q__1.real;
            i__2 = j;
            splus *= rhs[i__2].real;
            if(splus > sminu)
            {
                i__2 = j;
                rhs[i__2].real = bp.real;
                rhs[i__2].imag = bp.imag; // , expr subst
            }
            else if(sminu > splus)
            {
                i__2 = j;
                rhs[i__2].real = bm.real;
                rhs[i__2].imag = bm.imag; // , expr subst
            }
            else
            {
                /* In this case the updating sums are equal and we can */
                /* choose RHS(J) +1 or -1. The first time this happens we */
                /* choose -1, thereafter +1. This is a simple way to get */
                /* good estimates of matrices like Byers well-known example */
                /* (see [1]). (Not done in BSOLVE.) */
                i__2 = j;
                i__3 = j;
                q__1.real = rhs[i__3].real + pmone.real;
                q__1.imag = rhs[i__3].imag + pmone.imag; // , expr subst
                rhs[i__2].real = q__1.real;
                rhs[i__2].imag = q__1.imag; // , expr subst
                pmone.real = 1.f;
                pmone.imag = 0.f; // , expr subst
            }
            /* Compute the remaining r.h.s. */
            i__2 = j;
            q__1.real = -rhs[i__2].real;
            q__1.imag = -rhs[i__2].imag; // , expr subst
            temp.real = q__1.real;
            temp.imag = q__1.imag; // , expr subst
            i__2 = *n - j;
            aocl_blas_caxpy(&i__2, &temp, &z__[j + 1 + j * z_dim1], &c__1, &rhs[j + 1], &c__1);
            /* L10: */
        }
        /* Solve for U- part, lockahead for RHS(N) = +-1. This is not done */
        /* In BSOLVE and will hopefully give us a better estimate because */
        /* any ill-conditioning of the original matrix is transfered to U */
        /* and not to L. U(N, N) is an approximation to sigma_min(LU). */
        i__1 = *n - 1;
        aocl_blas_ccopy(&i__1, &rhs[1], &c__1, work, &c__1);
        i__1 = *n - 1;
        i__2 = *n;
        q__1.real = rhs[i__2].real + 1.f;
        q__1.imag = rhs[i__2].imag + 0.f; // , expr subst
        work[i__1].real = q__1.real;
        work[i__1].imag = q__1.imag; // , expr subst
        i__1 = *n;
        i__2 = *n;
        q__1.real = rhs[i__2].real - 1.f;
        q__1.imag = rhs[i__2].imag - 0.f; // , expr subst
        rhs[i__1].real = q__1.real;
        rhs[i__1].imag = q__1.imag; // , expr subst
        splus = 0.f;
        sminu = 0.f;
        for(i__ = *n; i__ >= 1; --i__)
        {
            c_div(&q__1, &c_b1, &z__[i__ + i__ * z_dim1]);
            temp.real = q__1.real;
            temp.imag = q__1.imag; // , expr subst
            i__1 = i__ - 1;
            i__2 = i__ - 1;
            q__1.real = work[i__2].real * temp.real - work[i__2].imag * temp.imag;
            q__1.imag = work[i__2].real * temp.imag + work[i__2].imag * temp.real; // , expr subst
            work[i__1].real = q__1.real;
            work[i__1].imag = q__1.imag; // , expr subst
            i__1 = i__;
            i__2 = i__;
            q__1.real = rhs[i__2].real * temp.real - rhs[i__2].imag * temp.imag;
            q__1.imag = rhs[i__2].real * temp.imag + rhs[i__2].imag * temp.real; // , expr subst
            rhs[i__1].real = q__1.real;
            rhs[i__1].imag = q__1.imag; // , expr subst
            i__1 = *n;
            for(k = i__ + 1; k <= i__1; ++k)
            {
                i__2 = i__ - 1;
                i__3 = i__ - 1;
                i__4 = k - 1;
                i__5 = i__ + k * z_dim1;
                q__3.real = z__[i__5].real * temp.real - z__[i__5].imag * temp.imag;
                q__3.imag = z__[i__5].real * temp.imag + z__[i__5].imag * temp.real; // , expr subst
                q__2.real = work[i__4].real * q__3.real - work[i__4].imag * q__3.imag;
                q__2.imag = work[i__4].real * q__3.imag + work[i__4].imag * q__3.real; // , expr subst
                q__1.real = work[i__3].real - q__2.real;
                q__1.imag = work[i__3].imag - q__2.imag; // , expr subst
                work[i__2].real = q__1.real;
                work[i__2].imag = q__1.imag; // , expr subst
                i__2 = i__;
                i__3 = i__;
                i__4 = k;
                i__5 = i__ + k * z_dim1;
                q__3.real = z__[i__5].real * temp.real - z__[i__5].imag * temp.imag;
                q__3.imag = z__[i__5].real * temp.imag + z__[i__5].imag * temp.real; // , expr subst
                q__2.real = rhs[i__4].real * q__3.real - rhs[i__4].imag * q__3.imag;
                q__2.imag = rhs[i__4].real * q__3.imag + rhs[i__4].imag * q__3.real; // , expr subst
                q__1.real = rhs[i__3].real - q__2.real;
                q__1.imag = rhs[i__3].imag - q__2.imag; // , expr subst
                rhs[i__2].real = q__1.real;
                rhs[i__2].imag = q__1.imag; // , expr subst
                /* L20: */
            }
            splus += c_abs(&work[i__ - 1]);
            sminu += c_abs(&rhs[i__]);
            /* L30: */
        }
        if(splus > sminu)
        {
            aocl_blas_ccopy(n, work, &c__1, &rhs[1], &c__1);
        }
        /* Apply the permutations JPIV to the computed solution (RHS) */
        i__1 = *n - 1;
        aocl_lapack_claswp(&c__1, &rhs[1], ldz, &c__1, &i__1, &jpiv[1], &c_n1);
        /* Compute the sum of squares */
        aocl_lapack_classq(n, &rhs[1], &c__1, rdscal, rdsum);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* ENTRY IJOB = 2 */
    /* Compute approximate nullvector XM of Z */
    aocl_lapack_cgecon("I", n, &z__[z_offset], ldz, &c_b24, &rtemp, work, rwork, &info);
    aocl_blas_ccopy(n, &work[*n], &c__1, xm, &c__1);
    /* Compute RHS */
    i__1 = *n - 1;
    aocl_lapack_claswp(&c__1, xm, ldz, &c__1, &i__1, &ipiv[1], &c_n1);
    aocl_lapack_cdotc_f2c(&q__3, n, xm, &c__1, xm, &c__1);
    c_sqrt(&q__2, &q__3);
    c_div(&q__1, &c_b1, &q__2);
    temp.real = q__1.real;
    temp.imag = q__1.imag; // , expr subst
    aocl_blas_cscal(n, &temp, xm, &c__1);
    aocl_blas_ccopy(n, xm, &c__1, xp, &c__1);
    aocl_blas_caxpy(n, &c_b1, &rhs[1], &c__1, xp, &c__1);
    q__1.real = -1.f;
    q__1.imag = -0.f; // , expr subst
    aocl_blas_caxpy(n, &q__1, xm, &c__1, &rhs[1], &c__1);
    aocl_lapack_cgesc2(n, &z__[z_offset], ldz, &rhs[1], &ipiv[1], &jpiv[1], &scale);
    aocl_lapack_cgesc2(n, &z__[z_offset], ldz, xp, &ipiv[1], &jpiv[1], &scale);
    if(aocl_blas_scasum(n, xp, &c__1) > aocl_blas_scasum(n, &rhs[1], &c__1))
    {
        aocl_blas_ccopy(n, xp, &c__1, &rhs[1], &c__1);
    }
    /* Compute the sum of squares */
    aocl_lapack_classq(n, &rhs[1], &c__1, rdscal, rdsum);
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLATDF */
}
/* clatdf_ */
