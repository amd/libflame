/* ../netlib/zlanhe.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__1 = 1;
 /* > \brief \b ZLANHE returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele ment of largest absolute value of a complex Hermitian matrix. */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download ZLANHE + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlanhe. f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlanhe. f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlanhe. f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* DOUBLE PRECISION FUNCTION ZLANHE( NORM, UPLO, N, A, LDA, WORK ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER NORM, UPLO */
 /* INTEGER LDA, N */
 /* .. */
 /* .. Array Arguments .. */
 /* DOUBLE PRECISION WORK( * ) */
 /* COMPLEX*16 A( LDA, * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > ZLANHE returns the value of the one norm, or the Frobenius norm, or */
 /* > the infinity norm, or the element of largest absolute value of a */
 /* > complex hermitian matrix A. */
 /* > \endverbatim */
 /* > */
 /* > \return ZLANHE */
 /* > \verbatim */
 /* > */
 /* > ZLANHE = ( max(f2c_abs(A(i,j))), NORM = 'M' or 'm' */
 /* > ( */
 /* > ( norm1(A), NORM = '1', 'O' or 'o' */
 /* > ( */
 /* > ( normI(A), NORM = 'I' or 'i' */
 /* > ( */
 /* > ( normF(A), NORM = 'F', 'f', 'E' or 'e' */
 /* > */
 /* > where norm1 denotes the one norm of a matrix (maximum column sum), */
 /* > normI denotes the infinity norm of a matrix (maximum row sum) and */
 /* > normF denotes the Frobenius norm of a matrix (square root of sum of */
 /* > squares). Note that max(f2c_abs(A(i,j))) is not a consistent matrix norm. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] NORM */
 /* > \verbatim */
 /* > NORM is CHARACTER*1 */
 /* > Specifies the value to be returned in ZLANHE as described */
 /* > above. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] UPLO */
 /* > \verbatim */
 /* > UPLO is CHARACTER*1 */
 /* > Specifies whether the upper or lower triangular part of the */
 /* > hermitian matrix A is to be referenced. */
 /* > = 'U': Upper triangular part of A is referenced */
 /* > = 'L': Lower triangular part of A is referenced */
 /* > \endverbatim */
 /* > */
 /* > \param[in] N */
 /* > \verbatim */
 /* > N is INTEGER */
 /* > The order of the matrix A. N >= 0. When N = 0, ZLANHE is */
 /* > set to zero. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] A */
 /* > \verbatim */
 /* > A is COMPLEX*16 array, dimension (LDA,N) */
 /* > The hermitian matrix A. If UPLO = 'U', the leading n by n */
 /* > upper triangular part of A contains the upper triangular part */
 /* > of the matrix A, and the strictly lower triangular part of A */
 /* > is not referenced. If UPLO = 'L', the leading n by n lower */
 /* > triangular part of A contains the lower triangular part of */
 /* > the matrix A, and the strictly upper triangular part of A is */
 /* > not referenced. Note that the imaginary parts of the diagonal */
 /* > elements need not be set and are assumed to be zero. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max(N,1). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)), */
 /* > where LWORK >= N when NORM = 'I' or '1' or 'O';
 otherwise, */
 /* > WORK is not referenced. */
 /* > \endverbatim */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date December 2016 */
 /* > \ingroup complex16HEauxiliary */
 /* ===================================================================== */
 doublereal zlanhe_(char *norm, char *uplo, integer *n, doublecomplex *a, integer *lda, doublereal *work) {
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE 
 char buffer[256]; 
 snprintf(buffer, 256,"zlanhe inputs: norm %c, uplo %c, n %d, lda %d",*norm, *uplo, *n, *lda);
 AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
 /* System generated locals */
 integer a_dim1, a_offset, i__1, i__2;
 doublereal ret_val, d__1;
 /* Builtin functions */
 double z_abs(doublecomplex *), sqrt(doublereal);
 /* Local variables */
 extern /* Subroutine */
 int dcombssq_(doublereal *, doublereal *);
 integer i__, j;
 doublereal sum, ssq[2], absa;
 extern logical lsame_(char *, char *);
 doublereal value;
 extern logical disnan_(doublereal *);
 doublereal colssq[2];
 extern /* Subroutine */
 int zlassq_(integer *, doublecomplex *, integer *, doublereal *, doublereal *);
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
 /* .. Local Arrays .. */
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
 --work;
 /* Function Body */
 if (*n == 0) {
 value = 0.;
 }
 else if (lsame_(norm, "M")) {
 /* Find max(f2c_abs(A(i,j))). */
 value = 0.;
 if (lsame_(uplo, "U")) {
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 i__2 = j - 1;
 for (i__ = 1;
 i__ <= i__2;
 ++i__) {
 sum = z_abs(&a[i__ + j * a_dim1]);
 if (value < sum || disnan_(&sum)) {
 value = sum;
 }
 /* L10: */
 }
 i__2 = j + j * a_dim1;
 sum = (d__1 = a[i__2].r, f2c_abs(d__1));
 if (value < sum || disnan_(&sum)) {
 value = sum;
 }
 /* L20: */
 }
 }
 else {
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 i__2 = j + j * a_dim1;
 sum = (d__1 = a[i__2].r, f2c_abs(d__1));
 if (value < sum || disnan_(&sum)) {
 value = sum;
 }
 i__2 = *n;
 for (i__ = j + 1;
 i__ <= i__2;
 ++i__) {
 sum = z_abs(&a[i__ + j * a_dim1]);
 if (value < sum || disnan_(&sum)) {
 value = sum;
 }
 /* L30: */
 }
 /* L40: */
 }
 }
 }
 else if (lsame_(norm, "I") || lsame_(norm, "O") || *(unsigned char *)norm == '1') {
 /* Find normI(A) ( = norm1(A), since A is hermitian). */
 value = 0.;
 if (lsame_(uplo, "U")) {
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 sum = 0.;
 i__2 = j - 1;
 for (i__ = 1;
 i__ <= i__2;
 ++i__) {
 absa = z_abs(&a[i__ + j * a_dim1]);
 sum += absa;
 work[i__] += absa;
 /* L50: */
 }
 i__2 = j + j * a_dim1;
 work[j] = sum + (d__1 = a[i__2].r, f2c_abs(d__1));
 /* L60: */
 }
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 sum = work[i__];
 if (value < sum || disnan_(&sum)) {
 value = sum;
 }
 /* L70: */
 }
 }
 else {
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 work[i__] = 0.;
 /* L80: */
 }
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 i__2 = j + j * a_dim1;
 sum = work[j] + (d__1 = a[i__2].r, f2c_abs(d__1));
 i__2 = *n;
 for (i__ = j + 1;
 i__ <= i__2;
 ++i__) {
 absa = z_abs(&a[i__ + j * a_dim1]);
 sum += absa;
 work[i__] += absa;
 /* L90: */
 }
 if (value < sum || disnan_(&sum)) {
 value = sum;
 }
 /* L100: */
 }
 }
 }
 else if (lsame_(norm, "F") || lsame_(norm, "E")) {
 /* Find normF(A). */
 /* SSQ(1) is scale */
 /* SSQ(2) is sum-of-squares */
 /* For better accuracy, sum each column separately. */
 ssq[0] = 0.;
 ssq[1] = 1.;
 /* Sum off-diagonals */
 if (lsame_(uplo, "U")) {
 i__1 = *n;
 for (j = 2;
 j <= i__1;
 ++j) {
 colssq[0] = 0.;
 colssq[1] = 1.;
 i__2 = j - 1;
 zlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, colssq, &colssq[1]);
 dcombssq_(ssq, colssq);
 /* L110: */
 }
 }
 else {
 i__1 = *n - 1;
 for (j = 1;
 j <= i__1;
 ++j) {
 colssq[0] = 0.;
 colssq[1] = 1.;
 i__2 = *n - j;
 zlassq_(&i__2, &a[j + 1 + j * a_dim1], &c__1, colssq, &colssq[ 1]);
 dcombssq_(ssq, colssq);
 /* L120: */
 }
 }
 ssq[1] *= 2;
 /* Sum diagonal */
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 i__2 = i__ + i__ * a_dim1;
 if (a[i__2].r != 0.) {
 i__2 = i__ + i__ * a_dim1;
 absa = (d__1 = a[i__2].r, f2c_abs(d__1));
 if (ssq[0] < absa) {
 /* Computing 2nd power */
 d__1 = ssq[0] / absa;
 ssq[1] = ssq[1] * (d__1 * d__1) + 1.;
 ssq[0] = absa;
 }
 else {
 /* Computing 2nd power */
 d__1 = absa / ssq[0];
 ssq[1] += d__1 * d__1;
 }
 }
 /* L130: */
 }
 value = ssq[0] * sqrt(ssq[1]);
 }
 ret_val = value;
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return ret_val;
 /* End of ZLANHE */
 }
 /* zlanhe_ */
 
