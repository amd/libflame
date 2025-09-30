/* ../netlib/v3.9.0/zlasyf_aa.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static dcomplex c_b6 = {{-1.}, {-0.}};
static aocl_int64_t c__1 = 1;
static dcomplex c_b8 = {{1.}, {0.}};
static dcomplex c_b19 = {{0.}, {0.}};
/* > \brief \b ZLASYF_AA */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLASYF_AA + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlasyf_
 * aa.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlasyf_
 * aa.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlasyf_
 * aa.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLASYF_AA( UPLO, J1, M, NB, A, LDA, IPIV, */
/* H, LDH, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER J1, M, NB, LDA, LDH */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), H( LDH, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLATRF_AA factorizes a panel of a scomplex symmetric matrix A using */
/* > the Aasen's algorithm. The panel consists of a set of NB rows of A */
/* > when UPLO is U, or a set of NB columns when UPLO is L. */
/* > */
/* > In order to factorize the panel, the Aasen's algorithm requires the */
/* > last row, or column, of the previous panel. The first row, or column, */
/* > of A is set to be the first row, or column, of an identity matrix, */
/* > which is used to factorize the first panel. */
/* > */
/* > The resulting J-th row of U, or J-th column of L, is stored in the */
/* > (J-1)-th row, or column, of A (without the unit diagonals), while */
/* > the diagonal and subdiagonal of A are overwritten by those of T. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
 */
/* > = 'L': Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] J1 */
/* > \verbatim */
/* > J1 is INTEGER */
/* > The location of the first row, or column, of the panel */
/* > within the submatrix of A, passed to this routine, e.g., */
/* > when called by ZSYTRF_AA, for the first panel, J1 is 1, */
/* > while for the remaining panels, J1 is 2. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The dimension of the submatrix. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > The dimension of the panel to be facotorized. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,M) for */
/* > the first panel, while dimension (LDA,M+1) for the */
/* > remaining panels. */
/* > */
/* > On entry, A contains the last row, or column, of */
/* > the previous panel, and the trailing submatrix of A */
/* > to be factorized, except for the first panel, only */
/* > the panel is passed. */
/* > */
/* > On exit, the leading panel is factorized. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (M) */
/* > Details of the row and column interchanges, */
/* > the row and column k were interchanged with the row and */
/* > column IPIV(k). */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* > H is COMPLEX*16 workspace, dimension (LDH,NB). */
/* > */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* > LDH is INTEGER */
/* > The leading dimension of the workspace H. LDH >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 workspace, dimension (M). */
/* > \endverbatim */
/* > */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2017 */
/* > \ingroup complex16SYcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zlasyf_aa_(char *uplo, aocl_int_t *j1, aocl_int_t *m, aocl_int_t *nb, dcomplex *a,
                aocl_int_t *lda, aocl_int_t *ipiv, dcomplex *h__, aocl_int_t *ldh,
                dcomplex *work)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zlasyf_aa(uplo, j1, m, nb, a, lda, ipiv, h__, ldh, work);
#else
    aocl_int64_t j1_64 = *j1;
    aocl_int64_t m_64 = *m;
    aocl_int64_t nb_64 = *nb;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldh_64 = *ldh;

    aocl_lapack_zlasyf_aa(uplo, &j1_64, &m_64, &nb_64, a, &lda_64, ipiv, h__, &ldh_64, work);
#endif
}

void aocl_lapack_zlasyf_aa(char *uplo, aocl_int64_t *j1, aocl_int64_t *m, aocl_int64_t *nb,
                           dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                           dcomplex *h__, aocl_int64_t *ldh, dcomplex *work)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlasyf_aa inputs: uplo %c, j1 %" FLA_IS ", m %" FLA_IS ", nb %" FLA_IS
                      ", lda %" FLA_IS ", ldh %" FLA_IS "",
                      *uplo, *j1, *m, *nb, *lda, *ldh);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, h_dim1, h_offset, i__1, i__2;
    dcomplex z__1;
    /* Builtin functions */
    void z_div(dcomplex *, dcomplex *, dcomplex *);
    /* Local variables */
    aocl_int64_t j, k, i1, k1, i2, mj;
    dcomplex piv, alpha;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    /* -- LAPACK computational routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2017 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --work;
    /* Function Body */
    j = 1;
    /* K1 is the first column of the panel to be factorized */
    /* i.e., K1 is 2 for the first block column, and 1 for the rest of the blocks */
    k1 = 2 - *j1 + 1;
    if(lsame_(uplo, "U", 1, 1))
    {
        /* ..................................................... */
        /* Factorize A as U**T*D*U using the upper triangle of A */
        /* ..................................................... */
    L10:
        if(j > fla_min(*m, *nb))
        {
            goto L20;
        }
        /* K is the column to be factorized */
        /* when being called from ZSYTRF_AA, */
        /* > for the first block column, J1 is 1, hence J1+J-1 is J, */
        /* > for the rest of the columns, J1 is 2, and J1+J-1 is J+1, */
        k = *j1 + j - 1;
        if(j == *m)
        {
            /* Only need to compute T(J, J) */
            mj = 1;
        }
        else
        {
            mj = *m - j + 1;
        }
        /* H(J:M, J) := A(J, J:M) - H(J:M, 1:(J-1)) * L(J1:(J-1), J), */
        /* where H(J:M, J) has been initialized to be A(J, J:M) */
        if(k > 2)
        {
            /* K is the column to be factorized */
            /* > for the first block column, K is J, skipping the first two */
            /* columns */
            /* > for the rest of the columns, K is J+1, skipping only the */
            /* first column */
            i__1 = j - k1;
            aocl_blas_zgemv("No transpose", &mj, &i__1, &c_b6, &h__[j + k1 * h_dim1], ldh,
                            &a[j * a_dim1 + 1], &c__1, &c_b8, &h__[j + j * h_dim1], &c__1);
        }
        /* Copy H(i:M, i) into WORK */
        aocl_blas_zcopy(&mj, &h__[j + j * h_dim1], &c__1, &work[1], &c__1);
        if(j > k1)
        {
            /* Compute WORK := WORK - L(J-1, J:M) * T(J-1,J), */
            /* where A(J-1, J) stores T(J-1, J) and A(J-2, J:M) stores U(J-1, J:M) */
            i__1 = k - 1 + j * a_dim1;
            z__1.r = -a[i__1].r;
            z__1.i = -a[i__1].i; // , expr subst
            alpha.r = z__1.r;
            alpha.i = z__1.i; // , expr subst
            aocl_blas_zaxpy(&mj, &alpha, &a[k - 2 + j * a_dim1], lda, &work[1], &c__1);
        }
        /* Set A(J, J) = T(J, J) */
        i__1 = k + j * a_dim1;
        a[i__1].r = work[1].r;
        a[i__1].i = work[1].i; // , expr subst
        if(j < *m)
        {
            /* Compute WORK(2:M) = T(J, J) L(J, (J+1):M) */
            /* where A(J, J) stores T(J, J) and A(J-1, (J+1):M) stores U(J, (J+1):M) */
            if(k > 1)
            {
                i__1 = k + j * a_dim1;
                z__1.r = -a[i__1].r;
                z__1.i = -a[i__1].i; // , expr subst
                alpha.r = z__1.r;
                alpha.i = z__1.i; // , expr subst
                i__1 = *m - j;
                aocl_blas_zaxpy(&i__1, &alpha, &a[k - 1 + (j + 1) * a_dim1], lda, &work[2], &c__1);
            }
            /* Find fla_max(|WORK(2:M)|) */
            i__1 = *m - j;
            i2 = aocl_blas_izamax(&i__1, &work[2], &c__1) + 1;
            i__1 = i2;
            piv.r = work[i__1].r;
            piv.i = work[i__1].i; // , expr subst
            /* Apply symmetric pivot */
            if(i2 != 2 && (piv.r != 0. || piv.i != 0.))
            {
                /* Swap WORK(I1) and WORK(I2) */
                i1 = 2;
                i__1 = i2;
                i__2 = i1;
                work[i__1].r = work[i__2].r;
                work[i__1].i = work[i__2].i; // , expr subst
                i__1 = i1;
                work[i__1].r = piv.r;
                work[i__1].i = piv.i; // , expr subst
                /* Swap A(I1, I1+1:M) with A(I1+1:M, I2) */
                i1 = i1 + j - 1;
                i2 = i2 + j - 1;
                i__1 = i2 - i1 - 1;
                aocl_blas_zswap(&i__1, &a[*j1 + i1 - 1 + (i1 + 1) * a_dim1], lda,
                                &a[*j1 + i1 + i2 * a_dim1], &c__1);
                /* Swap A(I1, I2+1:M) with A(I2, I2+1:M) */
                if(i2 < *m)
                {
                    i__1 = *m - i2;
                    aocl_blas_zswap(&i__1, &a[*j1 + i1 - 1 + (i2 + 1) * a_dim1], lda,
                                    &a[*j1 + i2 - 1 + (i2 + 1) * a_dim1], lda);
                }
                /* Swap A(I1, I1) with A(I2,I2) */
                i__1 = i1 + *j1 - 1 + i1 * a_dim1;
                piv.r = a[i__1].r;
                piv.i = a[i__1].i; // , expr subst
                i__1 = *j1 + i1 - 1 + i1 * a_dim1;
                i__2 = *j1 + i2 - 1 + i2 * a_dim1;
                a[i__1].r = a[i__2].r;
                a[i__1].i = a[i__2].i; // , expr subst
                i__1 = *j1 + i2 - 1 + i2 * a_dim1;
                a[i__1].r = piv.r;
                a[i__1].i = piv.i; // , expr subst
                /* Swap H(I1, 1:J1) with H(I2, 1:J1) */
                i__1 = i1 - 1;
                aocl_blas_zswap(&i__1, &h__[i1 + h_dim1], ldh, &h__[i2 + h_dim1], ldh);
                ipiv[i1] = (aocl_int_t)(i2);
                if(i1 > k1 - 1)
                {
                    /* Swap L(1:I1-1, I1) with L(1:I1-1, I2), */
                    /* skipping the first column */
                    i__1 = i1 - k1 + 1;
                    aocl_blas_zswap(&i__1, &a[i1 * a_dim1 + 1], &c__1, &a[i2 * a_dim1 + 1], &c__1);
                }
            }
            else
            {
                ipiv[j + 1] = (aocl_int_t)(j + 1);
            }
            /* Set A(J, J+1) = T(J, J+1) */
            i__1 = k + (j + 1) * a_dim1;
            a[i__1].r = work[2].r;
            a[i__1].i = work[2].i; // , expr subst
            if(j < *nb)
            {
                /* Copy A(J+1:M, J+1) into H(J:M, J), */
                i__1 = *m - j;
                aocl_blas_zcopy(&i__1, &a[k + 1 + (j + 1) * a_dim1], lda,
                                &h__[j + 1 + (j + 1) * h_dim1], &c__1);
            }
            /* Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1), */
            /* where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1) */
            if(j < *m - 1)
            {
                i__1 = k + (j + 1) * a_dim1;
                if(a[i__1].r != 0. || a[i__1].i != 0.)
                {
                    z_div(&z__1, &c_b8, &a[k + (j + 1) * a_dim1]);
                    alpha.r = z__1.r;
                    alpha.i = z__1.i; // , expr subst
                    i__1 = *m - j - 1;
                    aocl_blas_zcopy(&i__1, &work[3], &c__1, &a[k + (j + 2) * a_dim1], lda);
                    i__1 = *m - j - 1;
                    aocl_blas_zscal(&i__1, &alpha, &a[k + (j + 2) * a_dim1], lda);
                }
                else
                {
                    i__1 = *m - j - 1;
                    aocl_lapack_zlaset("Full", &c__1, &i__1, &c_b19, &c_b19,
                                       &a[k + (j + 2) * a_dim1], lda);
                }
            }
        }
        ++j;
        goto L10;
    L20:;
    }
    else
    {
        /* ..................................................... */
        /* Factorize A as L*D*L**T using the lower triangle of A */
        /* ..................................................... */
    L30:
        if(j > fla_min(*m, *nb))
        {
            goto L40;
        }
        /* K is the column to be factorized */
        /* when being called from ZSYTRF_AA, */
        /* > for the first block column, J1 is 1, hence J1+J-1 is J, */
        /* > for the rest of the columns, J1 is 2, and J1+J-1 is J+1, */
        k = *j1 + j - 1;
        if(j == *m)
        {
            /* Only need to compute T(J, J) */
            mj = 1;
        }
        else
        {
            mj = *m - j + 1;
        }
        /* H(J:M, J) := A(J:M, J) - H(J:M, 1:(J-1)) * L(J, J1:(J-1))^T, */
        /* where H(J:M, J) has been initialized to be A(J:M, J) */
        if(k > 2)
        {
            /* K is the column to be factorized */
            /* > for the first block column, K is J, skipping the first two */
            /* columns */
            /* > for the rest of the columns, K is J+1, skipping only the */
            /* first column */
            i__1 = j - k1;
            aocl_blas_zgemv("No transpose", &mj, &i__1, &c_b6, &h__[j + k1 * h_dim1], ldh,
                            &a[j + a_dim1], lda, &c_b8, &h__[j + j * h_dim1], &c__1);
        }
        /* Copy H(J:M, J) into WORK */
        aocl_blas_zcopy(&mj, &h__[j + j * h_dim1], &c__1, &work[1], &c__1);
        if(j > k1)
        {
            /* Compute WORK := WORK - L(J:M, J-1) * T(J-1,J), */
            /* where A(J-1, J) = T(J-1, J) and A(J, J-2) = L(J, J-1) */
            i__1 = j + (k - 1) * a_dim1;
            z__1.r = -a[i__1].r;
            z__1.i = -a[i__1].i; // , expr subst
            alpha.r = z__1.r;
            alpha.i = z__1.i; // , expr subst
            aocl_blas_zaxpy(&mj, &alpha, &a[j + (k - 2) * a_dim1], &c__1, &work[1], &c__1);
        }
        /* Set A(J, J) = T(J, J) */
        i__1 = j + k * a_dim1;
        a[i__1].r = work[1].r;
        a[i__1].i = work[1].i; // , expr subst
        if(j < *m)
        {
            /* Compute WORK(2:M) = T(J, J) L((J+1):M, J) */
            /* where A(J, J) = T(J, J) and A((J+1):M, J-1) = L((J+1):M, J) */
            if(k > 1)
            {
                i__1 = j + k * a_dim1;
                z__1.r = -a[i__1].r;
                z__1.i = -a[i__1].i; // , expr subst
                alpha.r = z__1.r;
                alpha.i = z__1.i; // , expr subst
                i__1 = *m - j;
                aocl_blas_zaxpy(&i__1, &alpha, &a[j + 1 + (k - 1) * a_dim1], &c__1, &work[2],
                                &c__1);
            }
            /* Find fla_max(|WORK(2:M)|) */
            i__1 = *m - j;
            i2 = aocl_blas_izamax(&i__1, &work[2], &c__1) + 1;
            i__1 = i2;
            piv.r = work[i__1].r;
            piv.i = work[i__1].i; // , expr subst
            /* Apply symmetric pivot */
            if(i2 != 2 && (piv.r != 0. || piv.i != 0.))
            {
                /* Swap WORK(I1) and WORK(I2) */
                i1 = 2;
                i__1 = i2;
                i__2 = i1;
                work[i__1].r = work[i__2].r;
                work[i__1].i = work[i__2].i; // , expr subst
                i__1 = i1;
                work[i__1].r = piv.r;
                work[i__1].i = piv.i; // , expr subst
                /* Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
                i1 = i1 + j - 1;
                i2 = i2 + j - 1;
                i__1 = i2 - i1 - 1;
                aocl_blas_zswap(&i__1, &a[i1 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1,
                                &a[i2 + (*j1 + i1) * a_dim1], lda);
                /* Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
                if(i2 < *m)
                {
                    i__1 = *m - i2;
                    aocl_blas_zswap(&i__1, &a[i2 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1,
                                    &a[i2 + 1 + (*j1 + i2 - 1) * a_dim1], &c__1);
                }
                /* Swap A(I1, I1) with A(I2, I2) */
                i__1 = i1 + (*j1 + i1 - 1) * a_dim1;
                piv.r = a[i__1].r;
                piv.i = a[i__1].i; // , expr subst
                i__1 = i1 + (*j1 + i1 - 1) * a_dim1;
                i__2 = i2 + (*j1 + i2 - 1) * a_dim1;
                a[i__1].r = a[i__2].r;
                a[i__1].i = a[i__2].i; // , expr subst
                i__1 = i2 + (*j1 + i2 - 1) * a_dim1;
                a[i__1].r = piv.r;
                a[i__1].i = piv.i; // , expr subst
                /* Swap H(I1, I1:J1) with H(I2, I2:J1) */
                i__1 = i1 - 1;
                aocl_blas_zswap(&i__1, &h__[i1 + h_dim1], ldh, &h__[i2 + h_dim1], ldh);
                ipiv[i1] = (aocl_int_t)(i2);
                if(i1 > k1 - 1)
                {
                    /* Swap L(1:I1-1, I1) with L(1:I1-1, I2), */
                    /* skipping the first column */
                    i__1 = i1 - k1 + 1;
                    aocl_blas_zswap(&i__1, &a[i1 + a_dim1], lda, &a[i2 + a_dim1], lda);
                }
            }
            else
            {
                ipiv[j + 1] = (aocl_int_t)(j + 1);
            }
            /* Set A(J+1, J) = T(J+1, J) */
            i__1 = j + 1 + k * a_dim1;
            a[i__1].r = work[2].r;
            a[i__1].i = work[2].i; // , expr subst
            if(j < *nb)
            {
                /* Copy A(J+1:M, J+1) into H(J+1:M, J), */
                i__1 = *m - j;
                aocl_blas_zcopy(&i__1, &a[j + 1 + (k + 1) * a_dim1], &c__1,
                                &h__[j + 1 + (j + 1) * h_dim1], &c__1);
            }
            /* Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1), */
            /* where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1) */
            if(j < *m - 1)
            {
                i__1 = j + 1 + k * a_dim1;
                if(a[i__1].r != 0. || a[i__1].i != 0.)
                {
                    z_div(&z__1, &c_b8, &a[j + 1 + k * a_dim1]);
                    alpha.r = z__1.r;
                    alpha.i = z__1.i; // , expr subst
                    i__1 = *m - j - 1;
                    aocl_blas_zcopy(&i__1, &work[3], &c__1, &a[j + 2 + k * a_dim1], &c__1);
                    i__1 = *m - j - 1;
                    aocl_blas_zscal(&i__1, &alpha, &a[j + 2 + k * a_dim1], &c__1);
                }
                else
                {
                    i__1 = *m - j - 1;
                    aocl_lapack_zlaset("Full", &i__1, &c__1, &c_b19, &c_b19, &a[j + 2 + k * a_dim1],
                                       lda);
                }
            }
        }
        ++j;
        goto L30;
    L40:;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZLASYF_AA */
}
/* zlasyf_aa__ */
