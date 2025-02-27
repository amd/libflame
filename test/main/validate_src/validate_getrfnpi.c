/******************************************************************************
 * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_getrfnpi.c
 *  @brief Defines validate function of GETRFNPI() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

/* Validates bottom left submatrix of array A */
#define validate_getrfnpi_test_L_BL_block(x, nrm_prefix, mat_type, realtype)               \
    realtype norm_A_BL = 0, norm = 0;                                                      \
    copy_matrix(datatype, "full", m_BL, n_BL, (mat_type *)A + nfact, lda, L_BL, m_BL);     \
    if(imatrix == 'O')                                                                     \
    {                                                                                      \
        for(int i = 0; i < n_BL; i++)                                                      \
        {                                                                                  \
            /* To handle large size values nrm2 is used */                                 \
            mat_type *vector = (mat_type *)L_BL + i * m_BL;                                \
            norm_A_BL = fla_test_max(norm_A_BL, nrm_prefix##nrm2_(&m_BL, vector, &i_one)); \
        }                                                                                  \
    }                                                                                      \
    else                                                                                   \
    {                                                                                      \
        norm_A_BL = fla_lapack_##x##lange("F", &m_BL, &n_BL, L_BL, &m_BL, NULL);           \
    }                                                                                      \
    /* (L_BL * U_TL) = A_BL */                                                             \
    x##trsm_("R", "U", "N", "N", &m_BL, &n_BL, &x##_one, A_test, &lda, L_BL, &m_BL);       \
    /* Compare value result and calculated value for L_BL */                               \
    matrix_difference(datatype, m_BL, n_BL, L_BL, m_BL, (mat_type *)A_test + nfact, lda);  \
    /* Calculate norm */                                                                   \
    if(imatrix == 'O')                                                                     \
    {                                                                                      \
        for(int i = 0; i < n_BL; i++)                                                      \
        {                                                                                  \
            /* To handle large size values nrm2 is used */                                 \
            mat_type *vector = (mat_type *)L_BL + i * m_BL;                                \
            norm = fla_test_max(norm, nrm_prefix##nrm2_(&m_BL, vector, &i_one));           \
        }                                                                                  \
    }                                                                                      \
    else                                                                                   \
    {                                                                                      \
        norm = fla_lapack_##x##lange("F", &m_BL, &n_BL, L_BL, &m_BL, NULL);                \
    }                                                                                      \
    resid3 = (norm / norm_A_BL) / (n_BL * eps);

/* Validates top right submatrix of array A */
#define validate_getrfnpi_test_U_TR_block(x, nrm_prefix, mat_type, realtype)                    \
    realtype norm_A_TR = 0, norm = 0;                                                           \
    copy_matrix(datatype, "full", m_TR, n_TR, (mat_type *)A + nfact * lda, lda, U_TR, m_TR);    \
    if(imatrix == 'O')                                                                          \
    {                                                                                           \
        for(int i = 0; i < n_TR; i++)                                                           \
        {                                                                                       \
            /* To handle large size values nrm2 is used */                                      \
            mat_type *vector = (mat_type *)U_TR + i * m_TR;                                     \
            norm_A_TR = fla_test_max(norm_A_TR, nrm_prefix##nrm2_(&m_TR, vector, &i_one));      \
        }                                                                                       \
    }                                                                                           \
    else                                                                                        \
    {                                                                                           \
        norm_A_TR = fla_lapack_##x##lange("F", &m_TR, &n_TR, U_TR, &m_TR, NULL);                \
    }                                                                                           \
    /* (L_TL * U_TR) = A_TR */                                                                  \
    x##trsm_("L", "L", "N", "U", &m_TR, &n_TR, &x##_one, A_test, &lda, U_TR, &m_TR);            \
    /* Compare value result and calculated value for U_TR */                                    \
    matrix_difference(datatype, m_TR, n_TR, U_TR, m_TR, (mat_type *)A_test + nfact * lda, lda); \
    /* Calculate norm */                                                                        \
    if(imatrix == 'O')                                                                          \
    {                                                                                           \
        for(int i = 0; i < n_TR; i++)                                                           \
        {                                                                                       \
            /* To handle large size values nrm2 is used */                                      \
            mat_type *vector = (mat_type *)U_TR + i * m_TR;                                     \
            norm = fla_test_max(norm, nrm_prefix##nrm2_(&m_TR, vector, &i_one));                \
        }                                                                                       \
    }                                                                                           \
    else                                                                                        \
    {                                                                                           \
        norm = fla_lapack_##x##lange("F", &m_TR, &n_TR, U_TR, &m_TR, NULL);                     \
    }                                                                                           \
    resid4 = (norm / norm_A_TR) / (n_TR * eps);

/* Validates if the bottom right submatrix has been correctly updated by
   getrfnpi */
#define validate_getrfnpi_test_BR_block(x, nrm_prefix, mat_type, realtype)                        \
    realtype norm_A_BR = 0, norm = 0;                                                             \
    copy_matrix(datatype, "full", m_BR, n_BR, (mat_type *)A + nfact * lda + nfact, lda, A_BR,     \
                m_BR);                                                                            \
    if(imatrix == 'O')                                                                            \
    {                                                                                             \
        for(int i = 0; i < n_BR; i++)                                                             \
        {                                                                                         \
            /* To handle large size values nrm2 is used */                                        \
            mat_type *vector = (mat_type *)A_BR + i * m_BR;                                       \
            norm_A_BR = fla_test_max(norm_A_BR, nrm_prefix##nrm2_(&m_BR, vector, &i_one));        \
        }                                                                                         \
    }                                                                                             \
    else                                                                                          \
    {                                                                                             \
        norm_A_BR = fla_lapack_##x##lange("F", &m_BR, &n_BR, A_BR, &m_BR, NULL);                  \
    }                                                                                             \
    /* (L_BR * U_BR) = A_BR - L_BL * U_TR */                                                      \
    x##gemm_("N", "N", &m_BR, &n_BR, &nfact, &x##_n_one, (mat_type *)A_test + nfact, &lda,        \
             (mat_type *)A_test + nfact * lda, &lda, &x##_one, A_BR, &m_BR);                      \
    /* Compare value result and calculated value for A_BR */                                      \
    matrix_difference(datatype, m_BR, n_BR, A_BR, m_BR, (mat_type *)A_test + nfact * lda + nfact, \
                      lda);                                                                       \
    /* Calculate norm */                                                                          \
    if(imatrix == 'O')                                                                            \
    {                                                                                             \
        for(int i = 0; i < n_BR; i++)                                                             \
        {                                                                                         \
            /* To handle large size values nrm2 is used */                                        \
            mat_type *vector = (mat_type *)A_BR + i * m_BR;                                       \
            norm = fla_test_max(norm, nrm_prefix##nrm2_(&m_BR, vector, &i_one));                  \
        }                                                                                         \
    }                                                                                             \
    else                                                                                          \
    {                                                                                             \
        norm = fla_lapack_##x##lange("F", &m_BR, &n_BR, A_BR, &m_BR, NULL);                       \
    }                                                                                             \
    resid5 = (norm / norm_A_BR) / (n_BR * eps);

#define validate_getrfnpi_run_tests(x, nrm_prefix, mat_type, realtype)        \
    /* If the size of bottom left submatrix is non zero, then validate        \
       if the bottom left submatrix has been correctly updated by getrfnpi */ \
    if(m_BL > 0 && n_BL > 0)                                                  \
    {                                                                         \
        create_matrix(datatype, LAPACK_COL_MAJOR, m_BL, n_BL, &L_BL, m_BL);   \
        validate_getrfnpi_test_L_BL_block(x, nrm_prefix, mat_type, realtype); \
        free_matrix(L_BL);                                                    \
    }                                                                         \
    /* If the size of top right submatrix is non zero, then validate          \
       if the bottom left submatrix has been correctly updated by getrfnpi */ \
    if(m_TR > 0 && n_TR > 0)                                                  \
    {                                                                         \
        create_matrix(datatype, LAPACK_COL_MAJOR, m_TR, n_TR, &U_TR, m_TR);   \
        validate_getrfnpi_test_U_TR_block(x, nrm_prefix, mat_type, realtype); \
        free_matrix(U_TR);                                                    \
    }                                                                         \
    /* If the size of bottom right submatrix is non zero, then validate       \
       if the bottom left submatrix has been correctly updated by getrfnpi */ \
    if(m_BR > 0 && n_BR > 0)                                                  \
    {                                                                         \
        create_matrix(datatype, LAPACK_COL_MAJOR, m_BR, n_BR, &A_BR, m_BR);   \
        validate_getrfnpi_test_BR_block(x, nrm_prefix, mat_type, realtype);   \
        free_matrix(A_BR);                                                    \
    }

void validate_getrfnpi(char *tst_api, integer m_A, integer n_A, integer nfact, void *A,
                       void *A_test, integer lda, integer *IPIV, integer datatype,
                       double err_thresh, char imatrix)
{
    /* System generated locals */
    integer m_BR, n_BR, m_BL, n_BL, m_TR, n_TR;
    void *A_BR, *L_BL, *U_TR;
    double residual, resid1 = 0., resid2 = 0.;
    double resid3 = 0., resid4 = 0., resid5 = 0.;

    /* Early return conditions */
    if(m_A == 0 || n_A == 0 || nfact == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m_A, n_A, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m_A, n_A, err_thresh);

    /* Call validate getrf on factored block */
    validate_getrf_internal(nfact, nfact, A, A_test, lda, IPIV, datatype, imatrix, &resid1,
                            &resid2);

    m_BR = m_A - nfact;
    n_BR = n_A - nfact;
    m_BL = m_BR;
    n_BL = nfact;
    m_TR = nfact;
    n_TR = n_BR;

    switch(datatype)
    {
        case FLOAT:
        {
            float eps = fla_lapack_slamch("Epsilon");
            validate_getrfnpi_run_tests(s, s, real, real);
            break;
        }
        case DOUBLE:
        {
            double eps = fla_lapack_dlamch("Epsilon");
            validate_getrfnpi_run_tests(d, d, doublereal, doublereal);
            break;
        }
        case COMPLEX:
        {
            float eps = fla_lapack_slamch("Epsilon");
            validate_getrfnpi_run_tests(c, sc, scomplex, real);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double eps = fla_lapack_dlamch("Epsilon");
            validate_getrfnpi_run_tests(z, dz, dcomplex, doublereal);
            break;
        }
        default:
            resid3 = resid4 = resid5 = 0.;
            break;
    }

    residual = fla_test_max(resid1, resid2);
    residual = fla_test_max(resid3, residual);
    residual = fla_test_max(resid4, residual);
    residual = fla_test_max(resid5, residual);

    FLA_PRINT_TEST_STATUS(m_A, n_A, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
    FLA_PRINT_SUBTEST_STATUS(resid4, err_thresh, "04");
    FLA_PRINT_SUBTEST_STATUS(resid5, err_thresh, "05");
}
