/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_gejsv.c
 *  @brief Defines validate function of GEJSV() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

#define gejsv_scale_singular_values(realtype_prefix, realtype)              \
    do                                                                      \
    {                                                                       \
        if(((realtype *)stat)[0] != ((realtype *)stat)[1])                  \
        {                                                                   \
            realtype scale = ((realtype *)stat)[0] / ((realtype *)stat)[1]; \
            realtype_prefix##scal_(&m, &scale, S_scaled, &i_one);           \
        }                                                                   \
    } while(0)

#define gejsv_validate_num_singular_values(realtype)   \
    {                                                  \
        integer num_svds = 0, i;                       \
        for(i = 0; i < n; i++)                         \
        {                                              \
            if(((realtype *)S_scaled)[i] > 0.0)        \
            {                                          \
                num_svds++;                            \
            }                                          \
        }                                              \
        resid1 = num_svds == istat[1] ? 0.0 : DBL_MAX; \
    }

#define gejsv_validate_singular_values(realtype_prefix, realtype)                        \
    {                                                                                    \
        realtype norm_orig = realtype_prefix##nrm2_(&n, S, &i_one);                      \
        void *S_copy = NULL;                                                             \
        create_vector(get_realtype(datatype), &S_copy, n);                               \
        copy_vector(get_realtype(datatype), n, S_scaled, 1, S_copy, 1);                  \
        realtype_prefix##axpy_(&n, &realtype_prefix##_n_one, S, &i_one, S_copy, &i_one); \
        realtype norm = realtype_prefix##nrm2_(&n, S_copy, &i_one);                      \
        resid2 = norm / (eps * norm_orig * n);                                           \
        free_vector(S_copy);                                                             \
    }

#define gejsv_validate_decomposition(type_prefix, realtype)                                       \
    {                                                                                             \
        realtype norm_orig = type_prefix##lange_("1", &m, &n, A, &lda, NULL);                     \
        void *A_copy = NULL;                                                                      \
        void *U_copy = NULL;                                                                      \
        create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_copy, lda);                            \
        create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &U_copy, ldu);                            \
        copy_matrix(datatype, "full", m, n, U, ldu, U_copy, ldu);                                 \
        multiply_matrix_diag_vector(datatype, 'R', VECTOR_TYPE_REAL, m, n, U_copy, ldu, S_scaled, \
                                    1);                                                           \
        fla_invoke_gemm(datatype, "N", "C", &m, &n, &n, U_copy, &ldu, V, &ldv, A_copy, &lda);     \
        matrix_difference(datatype, m, n, A_copy, lda, A, lda);                                   \
        realtype norm = type_prefix##lange_("1", &m, &n, A_copy, &lda, NULL);                     \
        resid3 = norm / (eps * norm_orig * n);                                                    \
        free_matrix(A_copy);                                                                      \
        free_matrix(U_copy);                                                                      \
    }

#define gejsv_validate_eliminated_singular_values                                                 \
    {                                                                                             \
        resid6 = 0.;                                                                              \
        /* If smaller noise values are not elimited then the calculated rank will be equal to the \
         * size of svd */                                                                         \
        if(istat[1] >= n)                                                                         \
        {                                                                                         \
            resid6 = DBL_MAX;                                                                     \
        }                                                                                         \
    }

#define gejsv_run_validations(type_prefix, realtype_prefix, realtype)                      \
    {                                                                                      \
        realtype eps = realtype_prefix##lamch_("P");                                       \
        gejsv_validate_num_singular_values(realtype);                                      \
        if(validate_singular_values)                                                       \
        {                                                                                  \
            gejsv_validate_singular_values(realtype_prefix, realtype);                     \
        }                                                                                  \
        else                                                                               \
        {                                                                                  \
            gejsv_validate_decomposition(type_prefix, realtype);                           \
        }                                                                                  \
        if(same_char(jobu, 'U') || same_char(jobu, 'F'))                                   \
        {                                                                                  \
            resid4 = check_orthogonality(datatype, U, m, svd_len, ldu);                    \
        }                                                                                  \
        if(same_char(jobv, 'V') || same_char(jobv, 'J'))                                   \
        {                                                                                  \
            resid5 = check_orthogonality(datatype, V, n, n, ldv);                          \
        }                                                                                  \
        /* If joba is A then validate that smaller singular values have been eliminated */ \
        if(test_eliminated_svds)                                                           \
        {                                                                                  \
            gejsv_validate_eliminated_singular_values;                                     \
        }                                                                                  \
    }

void validate_gejsv(char *tst_api, char joba, char jobu, char jobv, char jobr, char jobt, char jobp,
                    integer m, integer n, void *A, integer lda, void *S, void *S_test, void *U,
                    integer ldu, void *V, integer ldv, void *stat, integer *istat,
                    integer test_eliminated_svds, integer datatype, double err_thresh, void *scal,
                    char imatrix, void *params)
{
    double residual, resid1 = 0., resid2 = 0., resid3 = 0., resid4 = 0., resid5 = 0., resid6 = 0.;
    void *S_scaled = NULL;
    integer validate_singular_values;
    integer svd_len = same_char(jobu, 'F') ? m : n; /* Early return conditions */
    if(m == 0 || n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m, n, err_thresh);

    /* validate singular values only if U or V is not
    calculated */
    validate_singular_values = (same_char(jobu, 'N') || same_char(jobu, 'W') || same_char(jobv, 'N')
                                || same_char(jobv, 'W'));

    create_vector(get_realtype(datatype), &S_scaled, n);
    copy_vector(get_realtype(datatype), n, S_test, 1, S_scaled, 1);

    /* update sigular value by stat(0)/stat(1) */
    if(get_realtype(datatype) == FLOAT)
        gejsv_scale_singular_values(s, float);
    else
        gejsv_scale_singular_values(d, double);

    switch(datatype)
    {
        case FLOAT:
            gejsv_run_validations(s, s, float);
            break;
        case DOUBLE:
            gejsv_run_validations(d, d, double);
            break;
        case COMPLEX:
            gejsv_run_validations(c, s, float);
            break;
        case DOUBLE_COMPLEX:
            gejsv_run_validations(z, d, double);
            break;
    }

    free_vector(S_scaled);

    residual = fla_test_max(resid1, resid2);
    residual = fla_test_max(residual, resid3);
    residual = fla_test_max(residual, resid4);
    residual = fla_test_max(residual, resid5);
    residual = fla_test_max(residual, resid6);
    FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
    FLA_PRINT_SUBTEST_STATUS(resid4, err_thresh, "04");
    FLA_PRINT_SUBTEST_STATUS(resid5, err_thresh, "05");
    FLA_PRINT_SUBTEST_STATUS(resid6, err_thresh, "06");
}
