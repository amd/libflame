/* claqz1.f -- translated by f2c (version 20160102). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief \b CLAQZ1 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAQZ1 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/CLAQZ1.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/CLAQZ1.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/CLAQZ1.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAQZ1( ILQ, ILZ, K, ISTARTM, ISTOPM, IHI, A, LDA, B, */
/* $ LDB, NQ, QSTART, Q, LDQ, NZ, ZSTART, Z, LDZ ) */
/* IMPLICIT NONE */
/* Arguments */
/* LOGICAL, INTENT( IN ) :: ILQ, ILZ */
/* INTEGER, INTENT( IN ) :: K, LDA, LDB, LDQ, LDZ, ISTARTM, ISTOPM, */
/* $ NQ, NZ, QSTART, ZSTART, IHI */
/* COMPLEX :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAQZ1 chases a 1x1 shift bulge in a matrix pencil down a single position */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > */
/* > \param[in] ILQ */
/* > \verbatim */
/* > ILQ is LOGICAL */
/* > Determines whether or not to update the matrix Q */
/* > \endverbatim */
/* > */
/* > \param[in] ILZ */
/* > \verbatim */
/* > ILZ is LOGICAL */
/* > Determines whether or not to update the matrix Z */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > Index indicating the position of the bulge. */
/* > On entry, the bulge is located in */
/* > (A(k+1,k),B(k+1,k)). */
/* > On exit, the bulge is located in */
/* > (A(k+2,k+1),B(k+2,k+1)). */
/* > \endverbatim */
/* > */
/* > \param[in] ISTARTM */
/* > \verbatim */
/* > ISTARTM is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] ISTOPM */
/* > \verbatim */
/* > ISTOPM is INTEGER */
/* > Updates to (A,B) are restricted to */
/* > (istartm:k+2,k:istopm). It is assumed */
/* > without checking that istartm <= k+1 and */
/* > k+2 <= istopm */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* > IHI is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[inout] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of A as declared in */
/* > the calling procedure. */
/* > \endverbatim */
/* > \param[inout] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of B as declared in */
/* > the calling procedure. */
/* > \endverbatim */
/* > */
/* > \param[in] NQ */
/* > \verbatim */
/* > NQ is INTEGER */
/* > The order of the matrix Q */
/* > \endverbatim */
/* > */
/* > \param[in] QSTART */
/* > \verbatim */
/* > QSTART is INTEGER */
/* > Start index of the matrix Q. Rotations are applied */
/* > To columns k+2-qStart:k+3-qStart of Q. */
/* > \endverbatim */
/* > \param[inout] Q */
/* > \verbatim */
/* > Q is COMPLEX array, dimension (LDQ,NQ) */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of Q as declared in */
/* > the calling procedure. */
/* > \endverbatim */
/* > */
/* > \param[in] NZ */
/* > \verbatim */
/* > NZ is INTEGER */
/* > The order of the matrix Z */
/* > \endverbatim */
/* > */
/* > \param[in] ZSTART */
/* > \verbatim */
/* > ZSTART is INTEGER */
/* > Start index of the matrix Z. Rotations are applied */
/* > To columns k+1-qStart:k+2-qStart of Z. */
/* > \endverbatim */
/* > \param[inout] Z */
/* > \verbatim */
/* > Z is COMPLEX array, dimension (LDZ,NZ) */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of Q as declared in */
/* > the calling procedure. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Thijs Steel, KU Leuven */
/* > \date May 2020 */
/* > \ingroup complexGEcomputational */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void claqz1_(logical *ilq, logical *ilz, aocl_int_t *k, aocl_int_t *istartm, aocl_int_t *istopm,
             aocl_int_t *ihi, scomplex *a, aocl_int_t *lda, scomplex *b, aocl_int_t *ldb,
             aocl_int_t *nq, aocl_int_t *qstart, scomplex *q, aocl_int_t *ldq, aocl_int_t *nz,
             aocl_int_t *zstart, scomplex *z__, aocl_int_t *ldz)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_claqz1(ilq, ilz, k, istartm, istopm, ihi, a, lda, b, ldb, nq, qstart, q, ldq, nz,
                       zstart, z__, ldz);
#else
    aocl_int64_t k_64 = *k;
    aocl_int64_t istartm_64 = *istartm;
    aocl_int64_t istopm_64 = *istopm;
    aocl_int64_t ihi_64 = *ihi;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t nq_64 = *nq;
    aocl_int64_t qstart_64 = *qstart;
    aocl_int64_t ldq_64 = *ldq;
    aocl_int64_t nz_64 = *nz;
    aocl_int64_t zstart_64 = *zstart;
    aocl_int64_t ldz_64 = *ldz;

    aocl_lapack_claqz1(ilq, ilz, &k_64, &istartm_64, &istopm_64, &ihi_64, a, &lda_64, b, &ldb_64,
                       &nq_64, &qstart_64, q, &ldq_64, &nz_64, &zstart_64, z__, &ldz_64);
#endif
}

void aocl_lapack_claqz1(logical *ilq, logical *ilz, aocl_int64_t *k, aocl_int64_t *istartm,
                        aocl_int64_t *istopm, aocl_int64_t *ihi, scomplex *a, aocl_int64_t *lda,
                        scomplex *b, aocl_int64_t *ldb, aocl_int64_t *nq, aocl_int64_t *qstart,
                        scomplex *q, aocl_int64_t *ldq, aocl_int64_t *nz, aocl_int64_t *zstart,
                        scomplex *z__, aocl_int64_t *ldz)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, z_offset, i__1;
    scomplex q__1;
    /* Builtin functions */
    void r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    real c__;
    scomplex s, temp;
    extern void clartg_(scomplex *, scomplex *, real *, scomplex *, scomplex *);
    /* Arguments */
    /* Parameters */
    /* Local variables */
    /* External Functions */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    /* Function Body */
    if(*k + 1 == *ihi)
    {
        /* Shift is located on the edge of the matrix, remove it */
        clartg_(&b[*ihi + *ihi * b_dim1], &b[*ihi + (*ihi - 1) * b_dim1], &c__, &s, &temp);
        i__1 = *ihi + *ihi * b_dim1;
        b[i__1].real = temp.real;
        b[i__1].imag = temp.imag; // , expr subst
        i__1 = *ihi + (*ihi - 1) * b_dim1;
        b[i__1].real = 0.f;
        b[i__1].imag = 0.f; // , expr subst
        i__1 = *ihi - *istartm;
        aocl_lapack_crot(&i__1, &b[*istartm + *ihi * b_dim1], &c__1,
                         &b[*istartm + (*ihi - 1) * b_dim1], &c__1, &c__, &s);
        i__1 = *ihi - *istartm + 1;
        aocl_lapack_crot(&i__1, &a[*istartm + *ihi * a_dim1], &c__1,
                         &a[*istartm + (*ihi - 1) * a_dim1], &c__1, &c__, &s);
        if(*ilz)
        {
            aocl_lapack_crot(nz, &z__[(*ihi - *zstart + 1) * z_dim1 + 1], &c__1,
                             &z__[(*ihi - 1 - *zstart + 1) * z_dim1 + 1], &c__1, &c__, &s);
        }
    }
    else
    {
        /* Normal operation, move bulge down */
        /* Apply transformation from the right */
        clartg_(&b[*k + 1 + (*k + 1) * b_dim1], &b[*k + 1 + *k * b_dim1], &c__, &s, &temp);
        i__1 = *k + 1 + (*k + 1) * b_dim1;
        b[i__1].real = temp.real;
        b[i__1].imag = temp.imag; // , expr subst
        i__1 = *k + 1 + *k * b_dim1;
        b[i__1].real = 0.f;
        b[i__1].imag = 0.f; // , expr subst
        i__1 = *k + 2 - *istartm + 1;
        aocl_lapack_crot(&i__1, &a[*istartm + (*k + 1) * a_dim1], &c__1, &a[*istartm + *k * a_dim1],
                         &c__1, &c__, &s);
        i__1 = *k - *istartm + 1;
        aocl_lapack_crot(&i__1, &b[*istartm + (*k + 1) * b_dim1], &c__1, &b[*istartm + *k * b_dim1],
                         &c__1, &c__, &s);
        if(*ilz)
        {
            aocl_lapack_crot(nz, &z__[(*k + 1 - *zstart + 1) * z_dim1 + 1], &c__1,
                             &z__[(*k - *zstart + 1) * z_dim1 + 1], &c__1, &c__, &s);
        }
        /* Apply transformation from the left */
        clartg_(&a[*k + 1 + *k * a_dim1], &a[*k + 2 + *k * a_dim1], &c__, &s, &temp);
        i__1 = *k + 1 + *k * a_dim1;
        a[i__1].real = temp.real;
        a[i__1].imag = temp.imag; // , expr subst
        i__1 = *k + 2 + *k * a_dim1;
        a[i__1].real = 0.f;
        a[i__1].imag = 0.f; // , expr subst
        i__1 = *istopm - *k;
        aocl_lapack_crot(&i__1, &a[*k + 1 + (*k + 1) * a_dim1], lda, &a[*k + 2 + (*k + 1) * a_dim1],
                         lda, &c__, &s);
        i__1 = *istopm - *k;
        aocl_lapack_crot(&i__1, &b[*k + 1 + (*k + 1) * b_dim1], ldb, &b[*k + 2 + (*k + 1) * b_dim1],
                         ldb, &c__, &s);
        if(*ilq)
        {
            r_cnjg(&q__1, &s);
            aocl_lapack_crot(nq, &q[(*k + 1 - *qstart + 1) * q_dim1 + 1], &c__1,
                             &q[(*k + 2 - *qstart + 1) * q_dim1 + 1], &c__1, &c__, &q__1);
        }
    }
    /* End of CLAQZ1 */
    return;
}
/* claqz1_ */
