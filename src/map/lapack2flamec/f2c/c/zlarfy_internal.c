#include "FLA_f2c.h" /* Table of constant values */
static dcomplex c_b1 = {1., 0.};
static dcomplex c_b2 = {0., 0.};
static aocl_int64_t c__1 = 1;

void aocl_lapack_zlarfy(char *uplo, aocl_int64_t *n, dcomplex *v, aocl_int64_t *incv,
                        dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc,
                        dcomplex *work)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlarfy inputs: uplo %c, n %" FLA_IS ", incv %" FLA_IS ", ldc %" FLA_IS "",
                      *uplo, *n, *incv, *ldc);
    /* System generated locals */
    aocl_int64_t c_dim1, c_offset;
    dcomplex z__1, z__2, z__3, z__4;
    /* Local variables */
    dcomplex alpha;
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
    if(tau->real == 0. && tau->imag == 0.)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Form w:= C * v */
    aocl_blas_zhemv(uplo, n, &c_b1, &c__[c_offset], ldc, &v[1], incv, &c_b2, &work[1], &c__1);
    z__3.real = -.5;
    z__3.imag = -0.; // , expr subst
    z__2.real = z__3.real * tau->real - z__3.imag * tau->imag;
    z__2.imag = z__3.real * tau->imag + z__3.imag * tau->real; // , expr subst
    aocl_lapack_zdotc_f2c(&z__4, n, &work[1], &c__1, &v[1], incv);
    z__1.real = z__2.real * z__4.real - z__2.imag * z__4.imag;
    z__1.imag = z__2.real * z__4.imag + z__2.imag * z__4.real; // , expr subst
    alpha.real = z__1.real;
    alpha.imag = z__1.imag; // , expr subst
    aocl_blas_zaxpy(n, &alpha, &v[1], incv, &work[1], &c__1);
    /* C := C - v * w' - w * v' */
    z__1.real = -tau->real;
    z__1.imag = -tau->imag; // , expr subst
    aocl_blas_zher2(uplo, n, &z__1, &v[1], incv, &work[1], &c__1, &c__[c_offset], ldc);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZLARFY */
}
/* zlarfy_ */
