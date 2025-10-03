#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {1.f, 0.f};
static scomplex c_b2 = {0.f, 0.f};
static aocl_int64_t c__1 = 1;

void aocl_lapack_clarfy(char *uplo, aocl_int64_t *n, scomplex *v, aocl_int64_t *incv, scomplex *tau,
                        scomplex *c__, aocl_int64_t *ldc, scomplex *work)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clarfy inputs: uplo %c, n %lld, incv %lld, ldc %lld", *uplo, *n, *incv,
             *ldc);
#else
    snprintf(buffer, 256, "clarfy inputs: uplo %c, n %d, incv %d, ldc %d", *uplo, *n, *incv, *ldc);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t c_dim1, c_offset;
    scomplex q__1, q__2, q__3, q__4;
    /* Local variables */
    scomplex alpha;
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
    if(tau->real == 0.f && tau->imag == 0.f)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Form w:= C * v */
    aocl_blas_chemv(uplo, n, &c_b1, &c__[c_offset], ldc, &v[1], incv, &c_b2, &work[1], &c__1);
    q__3.real = -.5f;
    q__3.imag = -0.f; // , expr subst
    q__2.real = q__3.real * tau->real - q__3.imag * tau->imag;
    q__2.imag = q__3.real * tau->imag + q__3.imag * tau->real; // , expr subst
    aocl_lapack_cdotc_f2c(&q__4, n, &work[1], &c__1, &v[1], incv);
    q__1.real = q__2.real * q__4.real - q__2.imag * q__4.imag;
    q__1.imag = q__2.real * q__4.imag + q__2.imag * q__4.real; // , expr subst
    alpha.real = q__1.real;
    alpha.imag = q__1.imag; // , expr subst
    aocl_blas_caxpy(n, &alpha, &v[1], incv, &work[1], &c__1);
    /* C := C - v * w' - w * v' */
    q__1.real = -tau->real;
    q__1.imag = -tau->imag; // , expr subst
    aocl_blas_cher2(uplo, n, &q__1, &v[1], incv, &work[1], &c__1, &c__[c_offset], ldc);
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLARFY */
}
/* clarfy_ */
