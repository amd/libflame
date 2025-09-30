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
  SPFFRT2 computes partial or incomplete LDL' factorization of
  symmetric matrix in packed storage format.
*/

extern int sspffrt2_check(float *ap, aocl_int64_t *n, aocl_int64_t *ncolm, float *work, float *work2);
extern int dspffrt2_check(double *ap, aocl_int64_t *n, aocl_int64_t *ncolm, double *work, double *work2);
extern int cspffrt2_check(scomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm, scomplex *work,
                          scomplex *work2);
extern int zspffrt2_check(dcomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm, dcomplex *work,
                          dcomplex *work2);
extern void sspffrt2_fla(float *ap, aocl_int64_t *n, aocl_int64_t *ncolm, float *work, float *work2);
extern void dspffrt2_fla(double *ap, aocl_int64_t *n, aocl_int64_t *ncolm, double *work, double *work2);
extern void cspffrt2_fla(scomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm, scomplex *work, scomplex *work2);
extern void zspffrt2_fla(dcomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm, dcomplex *work, dcomplex *work2);

/** Generated wrapper function */
void sspffrt2_(real *ap, aocl_int_t *n, aocl_int_t *ncolm, real *work, real *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sspffrt2(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm_64 = *ncolm;

    aocl_lapack_sspffrt2(ap, &n_64, &ncolm_64, work, work2);
#endif
}

/** Generated wrapper function */
void dspffrt2_(doublereal *ap, aocl_int_t *n, aocl_int_t *ncolm, doublereal *work, doublereal *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dspffrt2(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm_64 = *ncolm;

    aocl_lapack_dspffrt2(ap, &n_64, &ncolm_64, work, work2);
#endif
}

/** Generated wrapper function */
void cspffrt2_(scomplex *ap, aocl_int_t *n, aocl_int_t *ncolm, scomplex *work, scomplex *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cspffrt2(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm_64 = *ncolm;

    aocl_lapack_cspffrt2(ap, &n_64, &ncolm_64, work, work2);
#endif
}

/** Generated wrapper function */
void zspffrt2_(dcomplex *ap, aocl_int_t *n, aocl_int_t *ncolm, dcomplex *work, dcomplex *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zspffrt2(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm_64 = *ncolm;

    aocl_lapack_zspffrt2(ap, &n_64, &ncolm_64, work, work2);
#endif
}

#define LAPACK_spffrt2(prefix)                                                        \
    void aocl_lapack_##prefix##spffrt2(PREFIX2LAPACK_TYPEDEF(prefix) * buff_AP, aocl_int64_t * n,  \
                               aocl_int64_t * ncolm, PREFIX2LAPACK_TYPEDEF(prefix) * work, \
                               PREFIX2LAPACK_TYPEDEF(prefix) * work2)

#define LAPACK_spffrt2_body(prefix) prefix##spffrt2_fla(buff_AP, n, ncolm, work, work2);

LAPACK_spffrt2(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sspffrt2 inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1(sspffrt2_check(buff_AP, n, ncolm, work, work2), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_spffrt2_body(s)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_spffrt2(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dspffrt2 inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1(dspffrt2_check(buff_AP, n, ncolm, work, work2), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_spffrt2_body(d)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_spffrt2(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cspffrt2 inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1(cspffrt2_check(buff_AP, n, ncolm, work, work2), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_spffrt2_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_spffrt2(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zspffrt2 inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1(zspffrt2_check(buff_AP, n, ncolm, work, work2), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_spffrt2_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}

#endif
