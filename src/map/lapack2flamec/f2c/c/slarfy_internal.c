#include "FLA_f2c.h" /* Table of constant values */
static real c_b2 = 1.f;
static real c_b3 = 0.f;
static aocl_int64_t c__1 = 1;

void aocl_lapack_slarfy(char *uplo, aocl_int64_t *n, real *v, aocl_int64_t *incv, real *tau,
                        real *c__, aocl_int64_t *ldc, real *work)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slarfy inputs: uplo %c ,n %" FLA_IS ",incv %" FLA_IS ",ldc %" FLA_IS "",
                      *uplo, *n, *incv, *ldc);
    /* System generated locals */
    aocl_int64_t c_dim1, c_offset;
    real r__1;
    /* Local variables */
    real alpha;
    /* -- LAPACK test routine (version 3.7.0) -- */
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    if(*tau == 0.f)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Form w:= C * v */
    aocl_blas_ssymv(uplo, n, &c_b2, &c__[c_offset], ldc, &v[1], incv, &c_b3, &work[1], &c__1);
    alpha = *tau * -.5f * aocl_blas_sdot(n, &work[1], &c__1, &v[1], incv);
    aocl_blas_saxpy(n, &alpha, &v[1], incv, &work[1], &c__1);
    /* C := C - v * w' - w * v' */
    r__1 = -(*tau);
    aocl_blas_ssyr2(uplo, n, &r__1, &v[1], incv, &work[1], &c__1, &c__[c_offset], ldc);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SLARFY */
}
/* slarfy_ */
