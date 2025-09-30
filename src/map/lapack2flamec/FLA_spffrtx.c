/*
    Copyright (c) 2020 Advanced Micro Devices, Inc.  All rights reserved.
    Oct 24, 2020
*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_prototypes.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_util_defs.h"

/*
  SPFFRTX computes partial or incomplete LDL' factorization of
  symmetric matrix in packed storage format.
*/

extern void sspffrtx_fla(float *ap, aocl_int64_t *n, aocl_int64_t *ncolm, float *work, float *work2);
extern void dspffrtx_fla(double *ap, aocl_int64_t *n, aocl_int64_t *ncolm, double *work, double *work2);
extern void cspffrtx_fla(scomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm, scomplex *work, scomplex *work2);
extern void zspffrtx_fla(dcomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm, dcomplex *work, dcomplex *work2);
extern int sspffrtx_check(float *ap, aocl_int64_t *n, aocl_int64_t *ncolm, float *work, float *work2);
extern int dspffrtx_check(double *ap, aocl_int64_t *n, aocl_int64_t *ncolm, double *work, double *work2);
extern int cspffrtx_check(scomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm, scomplex *work,
                          scomplex *work2);
extern int zspffrtx_check(dcomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm, dcomplex *work,
                          dcomplex *work2);

/** Generated wrapper function */
void sspffrtx_(real *ap, aocl_int_t *n, aocl_int_t *ncolm, real *work, real *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sspffrtx(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm_64 = *ncolm;

    aocl_lapack_sspffrtx(ap, &n_64, &ncolm_64, work, work2);
#endif
}

/** Generated wrapper function */
void dspffrtx_(doublereal *ap, aocl_int_t *n, aocl_int_t *ncolm, doublereal *work, doublereal *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dspffrtx(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm_64 = *ncolm;

    aocl_lapack_dspffrtx(ap, &n_64, &ncolm_64, work, work2);
#endif
}

/** Generated wrapper function */
void cspffrtx_(scomplex *ap, aocl_int_t *n, aocl_int_t *ncolm, scomplex *work, scomplex *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cspffrtx(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm_64 = *ncolm;

    aocl_lapack_cspffrtx(ap, &n_64, &ncolm_64, work, work2);
#endif
}

/** Generated wrapper function */
void zspffrtx_(dcomplex *ap, aocl_int_t *n, aocl_int_t *ncolm, dcomplex *work, dcomplex *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zspffrtx(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm_64 = *ncolm;

    aocl_lapack_zspffrtx(ap, &n_64, &ncolm_64, work, work2);
#endif
}

#define LAPACK_spffrtx(prefix)                                                        \
    void aocl_lapack_##prefix##spffrtx(PREFIX2LAPACK_TYPEDEF(prefix) * buff_AP, aocl_int64_t * n,  \
                               aocl_int64_t * ncolm, PREFIX2LAPACK_TYPEDEF(prefix) * work, \
                               PREFIX2LAPACK_TYPEDEF(prefix) * work2)

#define LAPACK_spffrtx_body(prefix) prefix##spffrtx_fla(buff_AP, n, ncolm, work, work2);

LAPACK_spffrtx(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sspffrtx inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1(sspffrtx_check(buff_AP, n, ncolm, work, work2), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_spffrtx_body(s)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_spffrtx(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dspffrtx inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1(dspffrtx_check(buff_AP, n, ncolm, work, work2), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_spffrtx_body(d)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_spffrtx(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cspffrtx inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1(cspffrtx_check(buff_AP, n, ncolm, work, work2), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_spffrtx_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_spffrtx(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zspffrtx inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1(zspffrtx_check(buff_AP, n, ncolm, work, work2), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_spffrtx_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}

#endif
