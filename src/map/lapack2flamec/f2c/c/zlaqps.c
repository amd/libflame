/* ../netlib/zlaqps.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static dcomplex c_b1 = {0., 0.};
static dcomplex c_b2 = {1., 0.};
static aocl_int64_t c__1 = 1;
/* > \brief \b ZLAQPS computes a step of QR factorization with column pivoting of a real m-by-n
 * matrix A by us ing BLAS level 3. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAQPS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqps.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqps.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqps.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1, */
/* VN2, AUXV, F, LDF ) */
/* .. Scalar Arguments .. */
/* INTEGER KB, LDA, LDF, M, N, NB, OFFSET */
/* .. */
/* .. Array Arguments .. */
/* INTEGER JPVT( * ) */
/* DOUBLE PRECISION VN1( * ), VN2( * ) */
/* COMPLEX*16 A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAQPS computes a step of QR factorization with column pivoting */
/* > of a scomplex M-by-N matrix A by using Blas-3. It tries to factorize */
/* > NB columns from A starting from the row OFFSET+1, and updates all */
/* > of the matrix with Blas-3 xGEMM. */
/* > */
/* > In some cases, due to catastrophic cancellations, it cannot */
/* > factorize NB columns. Hence, the actual number of factorized */
/* > columns is returned in KB. */
/* > */
/* > Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized. */
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
/* > \param[in] OFFSET */
/* > \verbatim */
/* > OFFSET is INTEGER */
/* > The number of rows of A that have been factorized in */
/* > previous steps. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > The number of columns to factorize. */
/* > \endverbatim */
/* > */
/* > \param[out] KB */
/* > \verbatim */
/* > KB is INTEGER */
/* > The number of columns actually factorized. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, block A(OFFSET+1:M,1:KB) is the triangular */
/* > factor obtained and block A(1:OFFSET,1:N) has been */
/* > accordingly pivoted, but no factorized. */
/* > The rest of the matrix, block A(OFFSET+1:M,KB+1:N) has */
/* > been updated. */
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
/* > JPVT(I) = K <==> Column K of the full matrix A has been */
/* > permuted into position I in AP. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX*16 array, dimension (KB) */
/* > The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VN1 */
/* > \verbatim */
/* > VN1 is DOUBLE PRECISION array, dimension (N) */
/* > The vector with the partial column norms. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VN2 */
/* > \verbatim */
/* > VN2 is DOUBLE PRECISION array, dimension (N) */
/* > The vector with the exact column norms. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AUXV */
/* > \verbatim */
/* > AUXV is COMPLEX*16 array, dimension (NB) */
/* > Auxiliar vector. */
/* > \endverbatim */
/* > */
/* > \param[in,out] F */
/* > \verbatim */
/* > F is COMPLEX*16 array, dimension (LDF,NB) */
/* > Matrix F**H = L * Y**H * A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDF */
/* > \verbatim */
/* > LDF is INTEGER */
/* > The leading dimension of the array F. LDF >= fla_max(1,N). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain */
/* > X. Sun, Computer Science Dept., Duke University, USA */
/* > \n */
/* > Partial column norm updating strategy modified on April 2011 */
/* > Z. Drmac and Z. Bujanovic, Dept. of Mathematics, */
/* > University of Zagreb, Croatia. */
/* > \par References: */
/* ================ */
/* > */
/* > LAPACK Working Note 176 */
/* > \htmlonly */
/* > <a href="http://www.netlib.org/lapack/lawnspdf/lawn176.pdf">[PDF]</a> */
/* > \endhtmlonly */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zlaqps_(aocl_int_t *m, aocl_int_t *n, aocl_int_t *offset, aocl_int_t *nb, aocl_int_t *kb,
             dcomplex *a, aocl_int_t *lda, aocl_int_t *jpvt, dcomplex *tau,
             doublereal *vn1, doublereal *vn2, dcomplex *auxv, dcomplex *f,
             aocl_int_t *ldf)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zlaqps(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t offset_64 = *offset;
    aocl_int64_t nb_64 = *nb;
    aocl_int64_t kb_64 = *kb;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldf_64 = *ldf;

    aocl_lapack_zlaqps(&m_64, &n_64, &offset_64, &nb_64, &kb_64, a, &lda_64, jpvt, tau, vn1, vn2,
                       auxv, f, &ldf_64);

    *kb = (aocl_int_t)kb_64;
#endif
}

void aocl_lapack_zlaqps(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *offset, aocl_int64_t *nb,
                        aocl_int64_t *kb, dcomplex *a, aocl_int64_t *lda, aocl_int_t *jpvt,
                        dcomplex *tau, doublereal *vn1, doublereal *vn2, dcomplex *auxv,
                        dcomplex *f, aocl_int64_t *ldf)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlaqps inputs: m %" FLA_IS ", n %" FLA_IS ", offset %" FLA_IS ", nb %" FLA_IS
                      ", kb %" FLA_IS ", lda %" FLA_IS ", ldf %" FLA_IS "",
                      *m, *n, *offset, *nb, *kb, *lda, *ldf);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, f_dim1, f_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    dcomplex z__1;
    /* Builtin functions */
    double sqrt(doublereal);
    void d_cnjg(dcomplex *, dcomplex *);
    double z_abs(dcomplex *);
    integer i_dnnt(doublereal *);
    /* Local variables */
    aocl_int64_t j, k, rk;
    dcomplex akk;
    aocl_int64_t pvt;
    doublereal temp, temp2, tol3z;
    aocl_int64_t itemp;
    extern doublereal dlamch_(char *);
    aocl_int64_t lsticc;
    aocl_int64_t lastrk;
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --jpvt;
    --tau;
    --vn1;
    --vn2;
    --auxv;
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    /* Function Body */
    /* Computing MIN */
    i__1 = *m;
    i__2 = *n + *offset; // , expr subst
    lastrk = fla_min(i__1, i__2);
    lsticc = 0;
    k = 0;
    tol3z = sqrt(dlamch_("Epsilon"));
    /* Beginning of while loop. */
L10:
    if(k < *nb && lsticc == 0)
    {
        ++k;
        rk = *offset + k;
        /* Determine ith pivot column and swap if necessary */
        i__1 = *n - k + 1;
        pvt = k - 1 + aocl_blas_idamax(&i__1, &vn1[k], &c__1);
        if(pvt != k)
        {
            aocl_blas_zswap(m, &a[pvt * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
            i__1 = k - 1;
            aocl_blas_zswap(&i__1, &f[pvt + f_dim1], ldf, &f[k + f_dim1], ldf);
            itemp = jpvt[pvt];
            jpvt[pvt] = jpvt[k];
            jpvt[k] = (aocl_int_t)(itemp);
            vn1[pvt] = vn1[k];
            vn2[pvt] = vn2[k];
        }
        /* Apply previous Householder reflectors to column K: */
        /* A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)**H. */
        if(k > 1)
        {
            i__1 = k - 1;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = k + j * f_dim1;
                d_cnjg(&z__1, &f[k + j * f_dim1]);
                f[i__2].real = z__1.real;
                f[i__2].imag = z__1.imag; // , expr subst
                /* L20: */
            }
            i__1 = *m - rk + 1;
            i__2 = k - 1;
            z__1.real = -1.;
            z__1.imag = -0.; // , expr subst
            aocl_blas_zgemv("No transpose", &i__1, &i__2, &z__1, &a[rk + a_dim1], lda,
                            &f[k + f_dim1], ldf, &c_b2, &a[rk + k * a_dim1], &c__1);
            i__1 = k - 1;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = k + j * f_dim1;
                d_cnjg(&z__1, &f[k + j * f_dim1]);
                f[i__2].real = z__1.real;
                f[i__2].imag = z__1.imag; // , expr subst
                /* L30: */
            }
        }
        /* Generate elementary reflector H(k). */
        if(rk < *m)
        {
            i__1 = *m - rk + 1;
            aocl_lapack_zlarfg(&i__1, &a[rk + k * a_dim1], &a[rk + 1 + k * a_dim1], &c__1, &tau[k]);
        }
        else
        {
            aocl_lapack_zlarfg(&c__1, &a[rk + k * a_dim1], &a[rk + k * a_dim1], &c__1, &tau[k]);
        }
        i__1 = rk + k * a_dim1;
        akk.real = a[i__1].real;
        akk.imag = a[i__1].imag; // , expr subst
        i__1 = rk + k * a_dim1;
        a[i__1].real = 1.;
        a[i__1].imag = 0.; // , expr subst
        /* Compute Kth column of F: */
        /* Compute F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**H*A(RK:M,K). */
        if(k < *n)
        {
            i__1 = *m - rk + 1;
            i__2 = *n - k;
            aocl_blas_zgemv("Conjugate transpose", &i__1, &i__2, &tau[k], &a[rk + (k + 1) * a_dim1],
                            lda, &a[rk + k * a_dim1], &c__1, &c_b1, &f[k + 1 + k * f_dim1], &c__1);
        }
        /* Padding F(1:K,K) with zeros. */
        i__1 = k;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j + k * f_dim1;
            f[i__2].real = 0.;
            f[i__2].imag = 0.; // , expr subst
            /* L40: */
        }
        /* Incremental updating of F: */
        /* F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)**H */
        /* *A(RK:M,K). */
        if(k > 1)
        {
            i__1 = *m - rk + 1;
            i__2 = k - 1;
            i__3 = k;
            z__1.real = -tau[i__3].real;
            z__1.imag = -tau[i__3].imag; // , expr subst
            aocl_blas_zgemv("Conjugate transpose", &i__1, &i__2, &z__1, &a[rk + a_dim1], lda,
                            &a[rk + k * a_dim1], &c__1, &c_b1, &auxv[1], &c__1);
            i__1 = k - 1;
            aocl_blas_zgemv("No transpose", n, &i__1, &c_b2, &f[f_dim1 + 1], ldf, &auxv[1], &c__1,
                            &c_b2, &f[k * f_dim1 + 1], &c__1);
        }
        /* Update the current row of A: */
        /* A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)**H. */
        if(k < *n)
        {
            i__1 = *n - k;
            z__1.real = -1.;
            z__1.imag = -0.; // , expr subst
            aocl_blas_zgemm("No transpose", "Conjugate transpose", &c__1, &i__1, &k, &z__1,
                            &a[rk + a_dim1], lda, &f[k + 1 + f_dim1], ldf, &c_b2,
                            &a[rk + (k + 1) * a_dim1], lda);
        }
        /* Update partial column norms. */
        if(rk < lastrk)
        {
            i__1 = *n;
            for(j = k + 1; j <= i__1; ++j)
            {
                if(vn1[j] != 0.)
                {
                    /* NOTE: The following 4 lines follow from the analysis in */
                    /* Lapack Working Note 176. */
                    temp = z_abs(&a[rk + j * a_dim1]) / vn1[j];
                    /* Computing MAX */
                    d__1 = 0.;
                    d__2 = (temp + 1.) * (1. - temp); // , expr subst
                    temp = fla_max(d__1, d__2);
                    /* Computing 2nd power */
                    d__1 = vn1[j] / vn2[j];
                    temp2 = temp * (d__1 * d__1);
                    if(temp2 <= tol3z)
                    {
                        vn2[j] = (doublereal)lsticc;
                        lsticc = j;
                    }
                    else
                    {
                        vn1[j] *= sqrt(temp);
                    }
                }
                /* L50: */
            }
        }
        i__1 = rk + k * a_dim1;
        a[i__1].real = akk.real;
        a[i__1].imag = akk.imag; // , expr subst
        /* End of while loop. */
        goto L10;
    }
    *kb = k;
    rk = *offset + *kb;
    /* Apply the block reflector to the rest of the matrix: */
    /* A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) - */
    /* A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)**H. */
    /* Computing MIN */
    i__1 = *n;
    i__2 = *m - *offset; // , expr subst
    if(*kb < fla_min(i__1, i__2))
    {
        i__1 = *m - rk;
        i__2 = *n - *kb;
        z__1.real = -1.;
        z__1.imag = -0.; // , expr subst
        aocl_blas_zgemm("No transpose", "Conjugate transpose", &i__1, &i__2, kb, &z__1,
                        &a[rk + 1 + a_dim1], lda, &f[*kb + 1 + f_dim1], ldf, &c_b2,
                        &a[rk + 1 + (*kb + 1) * a_dim1], lda);
    }
    /* Recomputation of difficult columns. */
L60:
    if(lsticc > 0)
    {
        itemp = i_dnnt(&vn2[lsticc]);
        i__1 = *m - rk;
        vn1[lsticc] = aocl_blas_dznrm2(&i__1, &a[rk + 1 + lsticc * a_dim1], &c__1);
        /* NOTE: The computation of VN1( LSTICC ) relies on the fact that */
        /* SNRM2 does not fail on vectors with norm below the value of */
        /* SQRT(DLAMCH('S')) */
        vn2[lsticc] = vn1[lsticc];
        lsticc = itemp;
        goto L60;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZLAQPS */
}
/* zlaqps_ */
