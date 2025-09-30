/* ../netlib/ctgexc.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CTGEXC */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CTGEXC + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctgexc.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctgexc.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctgexc.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, */
/* LDZ, IFST, ILST, INFO ) */
/* .. Scalar Arguments .. */
/* LOGICAL WANTQ, WANTZ */
/* INTEGER IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/* $ Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTGEXC reorders the generalized Schur decomposition of a scomplex */
/* > matrix pair (A,B), using an unitary equivalence transformation */
/* > (A, B) := Q * (A, B) * Z**H, so that the diagonal block of (A, B) with */
/* > row index IFST is moved to row ILST. */
/* > */
/* > (A, B) must be in generalized Schur canonical form, that is, A and */
/* > B are both upper triangular. */
/* > */
/* > Optionally, the matrices Q and Z of generalized Schur vectors are */
/* > updated. */
/* > */
/* > Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H */
/* > Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] WANTQ */
/* > \verbatim */
/* > WANTQ is LOGICAL */
/* > .TRUE. : update the left transformation matrix Q;
 */
/* > .FALSE.: do not update Q. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* > WANTZ is LOGICAL */
/* > .TRUE. : update the right transformation matrix Z;
 */
/* > .FALSE.: do not update Z. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the upper triangular matrix A in the pair (A, B). */
/* > On exit, the updated matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,N) */
/* > On entry, the upper triangular matrix B in the pair (A, B). */
/* > On exit, the updated matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is COMPLEX array, dimension (LDZ,N) */
/* > On entry, if WANTQ = .TRUE., the unitary matrix Q. */
/* > On exit, the updated matrix Q. */
/* > If WANTQ = .FALSE., Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. LDQ >= 1;
 */
/* > If WANTQ = .TRUE., LDQ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is COMPLEX array, dimension (LDZ,N) */
/* > On entry, if WANTZ = .TRUE., the unitary matrix Z. */
/* > On exit, the updated matrix Z. */
/* > If WANTZ = .FALSE., Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= 1;
 */
/* > If WANTZ = .TRUE., LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] IFST */
/* > \verbatim */
/* > IFST is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in,out] ILST */
/* > \verbatim */
/* > ILST is INTEGER */
/* > Specify the reordering of the diagonal blocks of (A, B). */
/* > The block with row index IFST is moved to row ILST, by a */
/* > sequence of swapping between adjacent blocks. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > =0: Successful exit. */
/* > <0: if INFO = -i, the i-th argument had an illegal value. */
/* > =1: The transformed matrix pair (A, B) would be too far */
/* > from generalized Schur form;
the problem is ill- */
/* > conditioned. (A, B) may have been partially reordered, */
/* > and ILST points to the first row of the current */
/* > position of the block being moved. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexGEcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* > Umea University, S-901 87 Umea, Sweden. */
/* > \par References: */
/* ================ */
/* > */
/* > [1] B. Kagstrom;
A Direct Method for Reordering Eigenvalues in the */
/* > Generalized Real Schur Form of a Regular Matrix Pair (A, B), in */
/* > M.S. Moonen et al (eds), Linear Algebra for Large Scale and */
/* > Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218. */
/* > \n */
/* > [2] B. Kagstrom and P. Poromaa;
Computing Eigenspaces with Specified */
/* > Eigenvalues of a Regular Matrix Pair (A, B) and Condition */
/* > Estimation: Theory, Algorithms and Software, Report */
/* > UMINF - 94.04, Department of Computing Science, Umea University, */
/* > S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87. */
/* > To appear in Numerical Algorithms, 1996. */
/* > \n */
/* > [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software */
/* > for Solving the Generalized Sylvester Equation and Estimating the */
/* > Separation between Regular Matrix Pairs, Report UMINF - 93.23, */
/* > Department of Computing Science, Umea University, S-901 87 Umea, */
/* > Sweden, December 1993, Revised April 1994, Also as LAPACK working */
/* > Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1, */
/* > 1996. */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void ctgexc_(logical *wantq, logical *wantz, aocl_int_t *n, scomplex *a, aocl_int_t *lda, scomplex *b,
             aocl_int_t *ldb, scomplex *q, aocl_int_t *ldq, scomplex *z__, aocl_int_t *ldz,
             aocl_int_t *ifst, aocl_int_t *ilst, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_ctgexc(wantq, wantz, n, a, lda, b, ldb, q, ldq, z__, ldz, ifst, ilst, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t ldq_64 = *ldq;
    aocl_int64_t ldz_64 = *ldz;
    aocl_int64_t ifst_64 = *ifst;
    aocl_int64_t ilst_64 = *ilst;
    aocl_int64_t info_64 = *info;

    aocl_lapack_ctgexc(wantq, wantz, &n_64, a, &lda_64, b, &ldb_64, q, &ldq_64, z__, &ldz_64,
                       &ifst_64, &ilst_64, &info_64);

    *ilst = (aocl_int_t)ilst_64;
    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_ctgexc(logical *wantq, logical *wantz, aocl_int64_t *n, scomplex *a,
                        aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, scomplex *q,
                        aocl_int64_t *ldq, scomplex *z__, aocl_int64_t *ldz, aocl_int64_t *ifst,
                        aocl_int64_t *ilst, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "ctgexc inputs: n %lld, lda %lld, ldb %lld, ldq %lld, ldz %lld, ifst %lld, ilst %lld",
             *n, *lda, *ldb, *ldq, *ldz, *ifst, *ilst);
#else
    snprintf(buffer, 256, "ctgexc inputs: n %d, lda %d, ldb %d, ldq %d, ldz %d, ifst %d, ilst %d",
             *n, *lda, *ldb, *ldq, *ldz, *ifst, *ilst);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, z_offset, i__1;
    /* Local variables */
    aocl_int64_t here;
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Decode and test input arguments. */
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
    *info = 0;
    if(*n < 0)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -5;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -7;
    }
    else if(*ldq < 1 || *wantq && *ldq < fla_max(1, *n))
    {
        *info = -9;
    }
    else if(*ldz < 1 || *wantz && *ldz < fla_max(1, *n))
    {
        *info = -11;
    }
    else if(*ifst < 1 || *ifst > *n)
    {
        *info = -12;
    }
    else if(*ilst < 1 || *ilst > *n)
    {
        *info = -13;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CTGEXC", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n <= 1)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(*ifst == *ilst)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(*ifst < *ilst)
    {
        here = *ifst;
    L10: /* Swap with next one below */
        aocl_lapack_ctgex2(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq,
                           &z__[z_offset], ldz, &here, info);
        if(*info != 0)
        {
            *ilst = here;
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
        ++here;
        if(here < *ilst)
        {
            goto L10;
        }
        --here;
    }
    else
    {
        here = *ifst - 1;
    L20: /* Swap with next one above */
        aocl_lapack_ctgex2(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq,
                           &z__[z_offset], ldz, &here, info);
        if(*info != 0)
        {
            *ilst = here;
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
        --here;
        if(here >= *ilst)
        {
            goto L20;
        }
        ++here;
    }
    *ilst = here;
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CTGEXC */
}
/* ctgexc_ */
