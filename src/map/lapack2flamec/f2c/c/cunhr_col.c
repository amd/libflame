/* ../netlib/v3.9.0/cunhr_col.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 = {1.f, 0.f};
static integer c__1 = 1;
/* > \brief \b CUNHR_COL */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CUNHR_COL + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunhr_c
 * ol.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunhr_c
 * ol.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunhr_c
 * ol.f"> */
/* > [TXT]</a> */
/* > */
/* Definition: */
/* =========== */
/* SUBROUTINE CUNHR_COL( M, N, NB, A, LDA, T, LDT, D, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDT, M, N, NB */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), D( * ), T( LDT, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CUNHR_COL takes an M-by-N complex matrix Q_in with orthonormal columns */
/* > as input, stored in A, and performs Householder Reconstruction (HR), */
/* > i.e. reconstructs Householder vectors V(i) implicitly representing */
/* > another M-by-N matrix Q_out, with the property that Q_in = Q_out*S, */
/* > where S is an N-by-N diagonal matrix with diagonal entries */
/* > equal to +1 or -1. The Householder vectors (columns V(i) of V) are */
/* > stored in A on output, and the diagonal entries of S are stored in D. */
/* > Block reflectors are also returned in T */
/* > (same output format as CGEQRT). */
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
/* > The number of columns of the matrix A. M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > The column block size to be used in the reconstruction */
/* > of Householder column vector blocks in the array A and */
/* > corresponding block reflectors in the array T. NB >= 1. */
/* > (Note that if NB > N, then N is used instead of NB */
/* > as the column block size.) */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > */
/* > On entry: */
/* > */
/* > The array A contains an M-by-N orthonormal matrix Q_in, */
/* > i.e the columns of A are orthogonal unit vectors. */
/* > */
/* > On exit: */
/* > */
/* > The elements below the diagonal of A represent the unit */
/* > lower-trapezoidal matrix V of Householder column vectors */
/* > V(i). The unit diagonal entries of V are not stored */
/* > (same format as the output below the diagonal in A from */
/* > CGEQRT). The matrix T and the matrix V stored on output */
/* > in A implicitly define Q_out. */
/* > */
/* > The elements above the diagonal contain the factor U */
/* > of the "modified" LU-decomposition: */
/* > Q_in - ( S ) = V * U */
/* > ( 0 ) */
/* > where 0 is a (M-N)-by-(M-N) zero matrix. */
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
/* > T is COMPLEX array, */
/* > dimension (LDT, N) */
/* > */
/* > Let NOCB = Number_of_output_col_blocks */
/* > = CEIL(N/NB) */
/* > */
/* > On exit, T(1:NB, 1:N) contains NOCB upper-triangular */
/* > block reflectors used to define Q_out stored in compact */
/* > form as a sequence of upper-triangular NB-by-NB column */
/* > blocks (same format as the output T in CGEQRT). */
/* > The matrix T and the matrix V stored on output in A */
/* > implicitly define Q_out. NOTE: The lower triangles */
/* > below the upper-triangular blcoks will be filled with */
/* > zeros. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. */
/* > LDT >= fla_max(1,fla_min(NB,N)). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is COMPLEX array, dimension fla_min(M,N). */
/* > The elements can be only plus or minus one. */
/* > */
/* > D(i) is constructed as D(i) = -SIGN(Q_in_i(i,i)), where */
/* > 1 <= i <= fla_min(M,N), and Q_in_i is Q_in after performing */
/* > i-1 steps of “modified” Gaussian elimination. */
/* > See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* > */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The computed M-by-M unitary factor Q_out is defined implicitly as */
/* > a product of unitary matrices Q_out(i). Each Q_out(i) is stored in */
/* > the compact WY-representation format in the corresponding blocks of */
/* > matrices V (stored in A) and T. */
/* > */
/* > The M-by-N unit lower-trapezoidal matrix V stored in the M-by-N */
/* > matrix A contains the column vectors V(i) in NB-size column */
/* > blocks VB(j). For example, VB(1) contains the columns */
/* > V(1), V(2), ... V(NB). NOTE: The unit entries on */
/* > the diagonal of Y are not stored in A. */
/* > */
/* > The number of column blocks is */
/* > */
/* > NOCB = Number_of_output_col_blocks = CEIL(N/NB) */
/* > */
/* > where each block is of order NB except for the last block, which */
/* > is of order LAST_NB = N - (NOCB-1)*NB. */
/* > */
/* > For example, if M=6, N=5 and NB=2, the matrix V is */
/* > */
/* > */
/* > V = ( VB(1), VB(2), VB(3) ) = */
/* > */
/* > = ( 1 ) */
/* > ( v21 1 ) */
/* > ( v31 v32 1 ) */
/* > ( v41 v42 v43 1 ) */
/* > ( v51 v52 v53 v54 1 ) */
/* > ( v61 v62 v63 v54 v65 ) */
/* > */
/* > */
/* > For each of the column blocks VB(i), an upper-triangular block */
/* > reflector TB(i) is computed. These blocks are stored as */
/* > a sequence of upper-triangular column blocks in the NB-by-N */
/* > matrix T. The size of each TB(i) block is NB-by-NB, except */
/* > for the last block, whose size is LAST_NB-by-LAST_NB. */
/* > */
/* > For example, if M=6, N=5 and NB=2, the matrix T is */
/* > */
/* > T = ( TB(1), TB(2), TB(3) ) = */
/* > */
/* > = ( t11 t12 t13 t14 t15 ) */
/* > ( t22 t24 ) */
/* > */
/* > */
/* > The M-by-M factor Q_out is given as a product of NOCB */
/* > unitary M-by-M matrices Q_out(i). */
/* > */
/* > Q_out = Q_out(1) * Q_out(2) * ... * Q_out(NOCB), */
/* > */
/* > where each matrix Q_out(i) is given by the WY-representation */
/* > using corresponding blocks from the matrices V and T: */
/* > */
/* > Q_out(i) = I - VB(i) * TB(i) * (VB(i))**T, */
/* > */
/* > where I is the identity matrix. Here is the formula with matrix */
/* > dimensions: */
/* > */
/* > Q(i){
M-by-M}
= I{
M-by-M}
- */
/* > VB(i){
M-by-INB}
* TB(i){
INB-by-INB}
* (VB(i))**T {
INB-by-M}
, */
/* > */
/* > where INB = NB, except for the last block NOCB */
/* > for which INB=LAST_NB. */
/* > */
/* > ===== */
/* > NOTE: */
/* > ===== */
/* > */
/* > If Q_in is the result of doing a QR factorization */
/* > B = Q_in * R_in, then: */
/* > */
/* > B = (Q_out*S) * R_in = Q_out * (S * R_in) = O_out * R_out. */
/* > */
/* > So if one wants to interpret Q_out as the result */
/* > of the QR factorization of B, then corresponding R_out */
/* > should be obtained by R_out = S * R_in, i.e. some rows of R_in */
/* > should be multiplied by -1. */
/* > */
/* > For the details of the algorithm, see [1]. */
/* > */
/* > [1] "Reconstructing Householder vectors from tall-skinny QR", */
/* > G. Ballard, J. Demmel, L. Grigori, M. Jacquelin, H.D. Nguyen, */
/* > E. Solomonik, J. Parallel Distrib. Comput., */
/* > vol. 85, pp. 3-31, 2015. */
/* > \endverbatim */
/* > */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2019 */
/* > \ingroup complexOTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > */
/* > November 2019, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
void cunhr_col_(integer *m, integer *n, integer *nb, complex *a, integer *lda, complex *t,
                integer *ldt, complex *d__, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256,
             "cunhr_col inputs: m %" FLA_IS ", n %" FLA_IS ", nb %" FLA_IS ", lda %" FLA_IS
             ", ldt %" FLA_IS "",
             *m, *n, *nb, *lda, *ldt);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1;
    /* Local variables */
    integer nplusone, i__, j, jb, jnb;
    extern /* Subroutine */
        void
        claunhr_col_getrfnp_(integer *, integer *, complex *, integer *, complex *, integer *),
        cscal_(integer *, complex *, complex *, integer *);
    integer iinfo;
    extern /* Subroutine */
        void
        ccopy_(integer *, complex *, integer *, complex *, integer *),
        ctrsm_(char *, char *, char *, char *, integer *, integer *, complex *, complex *,
               integer *, complex *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    integer jbtemp1, jbtemp2;
    /* -- LAPACK computational routine (version 3.9.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2019 */
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
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --d__;
    /* Function Body */
    *info = 0;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0 || *n > *m)
    {
        *info = -2;
    }
    else if(*nb < 1)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -5;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = fla_min(*nb, *n); // , expr subst
        if(*ldt < fla_max(i__1, i__2))
        {
            *info = -7;
        }
    }
    /* Handle error in the input parameters. */
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CUNHR_COL", &i__1, (ftnlen)9);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(fla_min(*m, *n) == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* On input, the M-by-N matrix A contains the unitary */
    /* M-by-N matrix Q_in. */
    /* (1) Compute the unit lower-trapezoidal V (ones on the diagonal */
    /* are not stored) by performing the "modified" LU-decomposition. */
    /* Q_in - ( S ) = V * U = ( V1 ) * U, */
    /* ( 0 ) ( V2 ) */
    /* where 0 is an (M-N)-by-N zero matrix. */
    /* (1-1) Factor V1 and U. */
    claunhr_col_getrfnp_(n, n, &a[a_offset], lda, &d__[1], &iinfo);
    /* (1-2) Solve for V2. */
    if(*m > *n)
    {
        i__1 = *m - *n;
        ctrsm_("R", "U", "N", "N", &i__1, n, &c_b1, &a[a_offset], lda, &a[*n + 1 + a_dim1], lda);
    }
    /* (2) Reconstruct the block reflector T stored in T(1:NB, 1:N) */
    /* as a sequence of upper-triangular blocks with NB-size column */
    /* blocking. */
    /* Loop over the column blocks of size NB of the array A(1:M,1:N) */
    /* and the array T(1:NB,1:N), JB is the column index of a column */
    /* block, JNB is the column block size at each step JB. */
    nplusone = *n + 1;
    i__1 = *n;
    i__2 = *nb;
    for(jb = 1; i__2 < 0 ? jb >= i__1 : jb <= i__1; jb += i__2)
    {
        /* (2-0) Determine the column block size JNB. */
        /* Computing MIN */
        i__3 = nplusone - jb;
        jnb = fla_min(i__3, *nb);
        /* (2-1) Copy the upper-triangular part of the current JNB-by-JNB */
        /* diagonal block U(JB) (of the N-by-N matrix U) stored */
        /* in A(JB:JB+JNB-1,JB:JB+JNB-1) into the upper-triangular part */
        /* of the current JNB-by-JNB block T(1:JNB,JB:JB+JNB-1) */
        /* column-by-column, total JNB*(JNB+1)/2 elements. */
        jbtemp1 = jb - 1;
        i__3 = jb + jnb - 1;
        for(j = jb; j <= i__3; ++j)
        {
            i__4 = j - jbtemp1;
            ccopy_(&i__4, &a[jb + j * a_dim1], &c__1, &t[j * t_dim1 + 1], &c__1);
        }
        /* (2-2) Perform on the upper-triangular part of the current */
        /* JNB-by-JNB diagonal block U(JB) (of the N-by-N matrix U) stored */
        /* in T(1:JNB,JB:JB+JNB-1) the following operation in place: */
        /* (-1)*U(JB)*S(JB), i.e the result will be stored in the upper- */
        /* triangular part of T(1:JNB,JB:JB+JNB-1). This multiplication */
        /* of the JNB-by-JNB diagonal block U(JB) by the JNB-by-JNB */
        /* diagonal block S(JB) of the N-by-N sign matrix S from the */
        /* right means changing the sign of each J-th column of the block */
        /* U(JB) according to the sign of the diagonal element of the block */
        /* S(JB), i.e. S(J,J) that is stored in the array element D(J). */
        i__3 = jb + jnb - 1;
        for(j = jb; j <= i__3; ++j)
        {
            i__4 = j;
            if(d__[i__4].r == 1.f && d__[i__4].i == 0.f)
            {
                i__4 = j - jbtemp1;
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cscal_(&i__4, &q__1, &t[j * t_dim1 + 1], &c__1);
            }
        }
        /* (2-3) Perform the triangular solve for the current block */
        /* matrix X(JB): */
        /* X(JB) * (A(JB)**T) = B(JB), where: */
        /* A(JB)**T is a JNB-by-JNB unit upper-triangular */
        /* coefficient block, and A(JB)=V1(JB), which */
        /* is a JNB-by-JNB unit lower-triangular block */
        /* stored in A(JB:JB+JNB-1,JB:JB+JNB-1). */
        /* The N-by-N matrix V1 is the upper part */
        /* of the M-by-N lower-trapezoidal matrix V */
        /* stored in A(1:M,1:N);
         */
        /* B(JB) is a JNB-by-JNB upper-triangular right-hand */
        /* side block, B(JB) = (-1)*U(JB)*S(JB), and */
        /* B(JB) is stored in T(1:JNB,JB:JB+JNB-1);
         */
        /* X(JB) is a JNB-by-JNB upper-triangular solution */
        /* block, X(JB) is the upper-triangular block */
        /* reflector T(JB), and X(JB) is stored */
        /* in T(1:JNB,JB:JB+JNB-1). */
        /* In other words, we perform the triangular solve for the */
        /* upper-triangular block T(JB): */
        /* T(JB) * (V1(JB)**T) = (-1)*U(JB)*S(JB). */
        /* Even though the blocks X(JB) and B(JB) are upper- */
        /* triangular, the routine CTRSM will access all JNB**2 */
        /* elements of the square T(1:JNB,JB:JB+JNB-1). Therefore, */
        /* we need to set to zero the elements of the block */
        /* T(1:JNB,JB:JB+JNB-1) below the diagonal before the call */
        /* to CTRSM. */
        /* (2-3a) Set the elements to zero. */
        jbtemp2 = jb - 2;
        i__3 = jb + jnb - 2;
        for(j = jb; j <= i__3; ++j)
        {
            i__4 = *nb;
            for(i__ = j - jbtemp2; i__ <= i__4; ++i__)
            {
                i__5 = i__ + j * t_dim1;
                t[i__5].r = 0.f;
                t[i__5].i = 0.f; // , expr subst
            }
        }
        /* (2-3b) Perform the triangular solve. */
        ctrsm_("R", "L", "C", "U", &jnb, &jnb, &c_b1, &a[jb + jb * a_dim1], lda,
               &t[jb * t_dim1 + 1], ldt);
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CUNHR_COL */
}
/* cunhr_col__ */
