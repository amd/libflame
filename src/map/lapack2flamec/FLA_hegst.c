/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_prototypes.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_util_defs.h"

/** Generated wrapper function */
void ssygst_(aocl_int_t *itype, char *uplo, aocl_int_t *m, real *buff_A, aocl_int_t *ldim_A, real *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_ssygst(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_ssygst(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dsygst_(aocl_int_t *itype, char *uplo, aocl_int_t *m, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dsygst(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dsygst(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void chegst_(aocl_int_t *itype, char *uplo, aocl_int_t *m, scomplex *buff_A, aocl_int_t *ldim_A, scomplex *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_chegst(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_chegst(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zhegst_(aocl_int_t *itype, char *uplo, aocl_int_t *m, dcomplex *buff_A, aocl_int_t *ldim_A, dcomplex *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zhegst(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zhegst(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void ssygs2_(aocl_int_t *itype, char *uplo, aocl_int_t *m, real *buff_A, aocl_int_t *ldim_A, real *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_ssygs2(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_ssygs2(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dsygs2_(aocl_int_t *itype, char *uplo, aocl_int_t *m, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dsygs2(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dsygs2(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void chegs2_(aocl_int_t *itype, char *uplo, aocl_int_t *m, scomplex *buff_A, aocl_int_t *ldim_A, scomplex *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_chegs2(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_chegs2(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zhegs2_(aocl_int_t *itype, char *uplo, aocl_int_t *m, dcomplex *buff_A, aocl_int_t *ldim_A, dcomplex *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zhegs2(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zhegs2(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

extern void zhegst_fla(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                       dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
extern void chegst_fla(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *b,
                       aocl_int64_t *ldb, aocl_int64_t *info);
extern void chegs2_fla(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *b,
                       aocl_int64_t *ldb, aocl_int64_t *info);
extern void zhegs2_fla(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                       dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);

/*
  ZHEGST reduces a scomplex Hermitian-definite generalized
  eigenproblem to standard form.
*/

#define LAPACK_hegst(prefix, name)                                                         \
    void aocl_lapack_##prefix##name##gst(aocl_int64_t *itype, char *uplo, aocl_int64_t *m,                   \
                                 PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, aocl_int64_t * ldim_A, \
                                 PREFIX2LAPACK_TYPEDEF(prefix) * buff_B, aocl_int64_t * ldim_B, \
                                 aocl_int64_t * info)

#define LAPACK_hegst_body(prefix)                              \
    FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);     \
    FLA_Inv inv_fla;                                           \
    FLA_Uplo uplo_fla;                                         \
    FLA_Obj A, B;                                              \
    FLA_Error init_result;                                     \
                                                               \
    FLA_Init_safe(&init_result);                               \
                                                               \
    FLA_Param_map_netlib_to_flame_inv((int *)itype, &inv_fla); \
    FLA_Param_map_netlib_to_flame_uplo(uplo, &uplo_fla);       \
                                                               \
    FLA_Obj_create_without_buffer(datatype, *m, *m, &A);       \
    FLA_Obj_attach_buffer(buff_A, 1, *ldim_A, &A);             \
                                                               \
    FLA_Obj_create_without_buffer(datatype, *m, *m, &B);       \
    FLA_Obj_attach_buffer(buff_B, 1, *ldim_B, &B);             \
                                                               \
    FLA_Eig_gest(inv_fla, uplo_fla, A, B);                     \
                                                               \
    FLA_Obj_free_without_buffer(&A);                           \
    FLA_Obj_free_without_buffer(&B);                           \
                                                               \
    FLA_Finalize_safe(init_result);                            \
                                                               \
    *info = 0;

LAPACK_hegst(s, sy)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("hegst-ssygst inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS "",
                      *itype, *uplo, *m, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1(ssygst_check(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(s)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_hegst(d, sy)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("hegst-dsygst inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS "",
                      *itype, *uplo, *m, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1(dsygst_check(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(d)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_hegst(c, he)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("chegst inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS "",
                      *itype, *uplo, *m, *ldim_A, *ldim_B);
#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1(chegst_check(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
#else
    {
        chegst_fla(itype, uplo, m, (scomplex *)buff_A, ldim_A, (scomplex *)buff_B, ldim_B, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#endif
}
LAPACK_hegst(z, he)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhegst inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS "",
                      *itype, *uplo, *m, *ldim_A, *ldim_B);
#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1(zhegst_check(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
#else
    {
        zhegst_fla(itype, uplo, m, (dcomplex *)buff_A, ldim_A, (dcomplex *)buff_B, ldim_B,
                   info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#endif
}

#define LAPACK_hegs2(prefix, name)                                                         \
    void aocl_lapack_##prefix##name##gs2(aocl_int64_t *itype, char *uplo, aocl_int64_t *m,                   \
                                 PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, aocl_int64_t * ldim_A, \
                                 PREFIX2LAPACK_TYPEDEF(prefix) * buff_B, aocl_int64_t * ldim_B, \
                                 aocl_int64_t * info)

LAPACK_hegs2(s, sy)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("hegs2-ssygs2 inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS "",
                      *itype, *uplo, *m, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1(ssygs2_check(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(s)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_hegs2(d, sy)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("hegs2-dsygs2 inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS "",
                      *itype, *uplo, *m, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1(dsygs2_check(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(d)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_hegs2(c, he)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("chegs2 inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS "",
                      *itype, *uplo, *m, *ldim_A, *ldim_B);
#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1(chegs2_check(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
#else
    {
        chegs2_fla(itype, uplo, m, (scomplex *)buff_A, ldim_A, (scomplex *)buff_B, ldim_B, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#endif
}
LAPACK_hegs2(z, he)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhegs2 inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS "",
                      *itype, *uplo, *m, *ldim_A, *ldim_B);
#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1(zhegs2_check(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
#else
    {
        zhegs2_fla(itype, uplo, m, (dcomplex *)buff_A, ldim_A, (dcomplex *)buff_B, ldim_B,
                   info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#endif
}

#endif
