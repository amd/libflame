/*
    Copyright (C) 2024-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/* > \brief \b validate_syev.c                                              */
/* =========== DOCUMENTATION ===========                                     */
/* Definition:                                                               */
/* ===========                                                               */
/* SUBROUTINE  validate_syev(&jobz, &range, n, A, A_test, lda, vl, vu, il,   */
/*                           iu, L, w, datatype, residual);                  */
/* > \par Purpose:                                                           */
/* =============                                                             */
/* >                                                                         */
/* > \verbatim                                                               */
/* >                                                                         */
/* > Defines validate function of symmetric eigen APIs to use in test suite  */
/* >                                                                         */
/* > This API validates symmetric eigen APIs functionality by performing     */
/* > the folllowing tests:                                                   */
/* >                                                                         */
/* > TEST 1: Verify A = (Z * lambda * Z') using norm calculation             */
/* >                                                                         */
/* > TEST 2: Check orthogonality of eigen vectors (Z * Z' = I)               */
/* >                                                                         */
/* > Test 3: check if any of the eigen vectors failed to converge.           */
/* >         Note: When eigen vectors fail to converge the corresponding     */
/* >               element in ifail != 0                                     */
/* >                                                                         */
/* > TEST 4: verify the input and output EVs are same or not                 */
/* >                                                                         */
/* > TEST I: When range = I, verify if EVs between the speficied indices     */
/* >         il, iu are generated                                            */
/* >                                                                         */
/* > \endverbatim                                                            */
/* ========================================================================= */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_syev(char *tst_api, char *jobz, char *range, integer n, void *A, void *A_test,
                   integer lda, integer il, integer iu, void *L, void *lambda, void *ifail,
                   integer datatype, double err_thresh, char imatrix, void *scal)
{
    double residual, resid1 = 0., resid2 = 0.;
    double resid3 = 0., resid4 = 0., resid5 = 0.;

    /* Early return conditions */
    if(n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(n, n, err_thresh);

    if(L != NULL)
    {
        sort_realtype_vector(datatype, "A", n, L, 1);
    }

    if(same_char(*range, 'I'))
    {
        /* Test I range
           check if output EVs matches the input EVs in given index range */
        if((L != NULL) && compare_realtype_vector(datatype, (iu - il + 1), L, 1, il, lambda, 1))
        {
            resid1 = DBL_MAX;
        }
    }
    else /* range A or V */
    {
        if(!same_char(*jobz, 'N'))
        {
            void *Z = NULL, *work = NULL, *Q = NULL;
            integer i;
            integer *buff = (integer *)ifail;

            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z, lda);
            reset_matrix(datatype, n, n, Z, lda);
            copy_matrix(datatype, "full", n, n, A_test, lda, Z, lda);

            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q, lda);
            reset_matrix(datatype, n, n, Q, lda);
            copy_matrix(datatype, "full", n, n, A_test, lda, Q, lda);

            /* Multiply Q * lambda(eigen values) */
            multiply_matrix_diag_vector(datatype, n, n, Q, lda, lambda, 1);

            switch(datatype)
            {
                case FLOAT:
                {
                    float norm, norm_A, eps;
                    eps = fla_lapack_slamch("P");

                    /* Test 2
                       compute norm(A - ((Q * lambda) * Q')) /
                                    (n * norm(A) * EPS)      */
                    norm_A = fla_lapack_slange("1", &n, &n, A, &lda, work);
                    sgemm_("N", "T", &n, &n, &n, &s_one, Q, &lda, Z, &lda, &s_n_one, A, &lda);
                    norm = fla_lapack_slange("1", &n, &n, A, &lda, work);
                    resid2 = norm / (eps * norm_A * (float)n);

                    /* Test 3
                       compute norm(I - Z'*Z) / (N * EPS)  */
                    resid3 = (float)check_orthogonality(datatype, Z, n, n, lda);
                    break;
                }

                case DOUBLE:
                {
                    double norm, norm_A, eps;
                    eps = fla_lapack_dlamch("P");

                    /* Test 2
                       compute norm(A - (Q * lambda * Q')) /
                                     (n * norm(A) * EPS)       */
                    norm_A = fla_lapack_dlange("1", &n, &n, A, &lda, work);
                    dgemm_("N", "T", &n, &n, &n, &d_one, Q, &lda, Z, &lda, &d_n_one, A, &lda);
                    norm = fla_lapack_dlange("1", &n, &n, A, &lda, work);
                    resid2 = norm / (eps * norm_A * (double)n);

                    /* Test 3
                       compute norm(I - Z'*Z) / (N * EPS)*/
                    resid3 = check_orthogonality(datatype, Z, n, n, lda);
                    break;
                }

                case COMPLEX:
                {
                    float norm, norm_A, eps;
                    eps = fla_lapack_slamch("P");

                    /* Test 2
                       compute norm(A - (Q * lambda * Q')) /
                                       (n * norm(A) * EPS)   */
                    norm_A = fla_lapack_clange("1", &n, &n, A, &lda, work);
                    cgemm_("N", "C", &n, &n, &n, &c_one, Q, &lda, Z, &lda, &c_n_one, A, &lda);
                    norm = fla_lapack_clange("1", &n, &n, A, &lda, work);
                    resid2 = norm / (eps * norm_A * (float)n);

                    /* Test 3
                       compute norm(I - Z'*Z) / (N * EPS)*/
                    resid3 = (float)check_orthogonality(datatype, Z, n, n, lda);
                    break;
                }

                case DOUBLE_COMPLEX:
                {
                    double norm, norm_A, eps;
                    eps = fla_lapack_dlamch("P");

                    /* Test 2
                       compute norm(A - (Q * lambda * Q')) /
                                    (n * norm(A) * EPS)      */
                    norm_A = fla_lapack_zlange("1", &n, &n, A, &lda, work);
                    zgemm_("N", "C", &n, &n, &n, &z_one, Q, &lda, Z, &lda, &z_n_one, A, &lda);
                    norm = fla_lapack_zlange("1", &n, &n, A, &lda, work);
                    resid2 = norm / (eps * norm_A * (double)n);

                    /* Test 3
                       compute norm(I - Z'*Z) / (N * EPS)*/
                    resid3 = check_orthogonality(datatype, Z, n, n, lda);
                    break;
                }
            }
            /* Test 4
               check if any of the eigen vectors failed to converge.
               Note: When eigen vectors fail to converge the corresponding
               element in ifail != 0 */
            if(ifail != NULL)
            {
                for(i = 0; i < n; i++)
                {
                    if(buff[i] != 0)
                        resid4 = DBL_MAX;
                }
            }

            free_matrix(Q);
            free_matrix(Z);
        }
        /* Test 5: In case of specific input generation, compare input and
           output eigen values */
        if(L != NULL)
        {
            void *work = NULL;
            if(get_realtype(datatype) == FLOAT)
            {
                float norm, norm_L, eps;
                eps = fla_lapack_slamch("P");

                if((same_char(imatrix, 'O') || same_char(imatrix, 'U')) && (scal != NULL))
                {
                    sscal_(&n, scal, L, &i_one);
                }

                norm_L = fla_lapack_slange("1", &n, &i_one, L, &i_one, work);
                saxpy_(&n, &s_n_one, lambda, &i_one, L, &i_one);
                norm = fla_lapack_slange("1", &n, &i_one, L, &i_one, work);
                resid5 = norm / (eps * norm_L * n);
            }
            else
            {
                double norm, norm_L, eps;
                eps = fla_lapack_dlamch("P");

                if((same_char(imatrix, 'O') || same_char(imatrix, 'U')) && (scal != NULL))
                {
                    dscal_(&n, scal, L, &i_one);
                }
                norm_L = fla_lapack_dlange("1", &n, &i_one, L, &i_one, work);
                daxpy_(&n, &d_n_one, lambda, &i_one, L, &i_one);
                norm = fla_lapack_dlange("1", &n, &i_one, L, &i_one, work);
                resid5 = norm / (eps * norm_L * n);
            }
        }
    }

    residual = fla_test_max(resid1, resid2);
    residual = fla_test_max(resid3, residual);
    residual = fla_test_max(resid4, residual);
    residual = fla_test_max(resid5, residual);

    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
    FLA_PRINT_SUBTEST_STATUS(resid4, err_thresh, "04");
    FLA_PRINT_SUBTEST_STATUS(resid5, err_thresh, "05");
}
