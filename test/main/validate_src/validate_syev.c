/******************************************************************************
 * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

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

void validate_syev(char *jobz, char *range, integer n, void *A, void *A_test, integer lda,
                   integer il, integer iu, void *L, void *lambda, void *ifail, integer datatype,
                   double *residual, char imatrix, void *scal)
{
    if(n == 0)
        return;
    *residual = 0;

    if(L != NULL)
    {
        sort_realtype_vector(datatype, "A", n, L, 1);
    }

    if((*range == 'I') || (*range == 'i'))
    {
        /* Test I range
           check if output EVs matches the input EVs in given index range */
        if((L != NULL) && compare_realtype_vector(datatype, (iu - il + 1), L, 1, il, lambda, 1))
        {
            *residual = DBL_MAX;
        }
    }
    else /* range A or V */
    {
        if(*jobz != 'N')
        {
            void *Z = NULL, *work = NULL, *Q = NULL;
            integer i;
            integer *buff = (integer *)ifail;

            create_matrix(datatype, matrix_layout, n, n, &Z, lda);
            reset_matrix(datatype, n, n, Z, lda);
            copy_matrix(datatype, "full", n, n, A_test, lda, Z, lda);

            create_matrix(datatype, matrix_layout, n, n, &Q, lda);
            reset_matrix(datatype, n, n, Q, lda);
            copy_matrix(datatype, "full", n, n, A_test, lda, Q, lda);

            /* Multiply Q * lambda(eigen values) */
            multiply_matrix_diag_vector(datatype, n, n, Q, lda, lambda, 1);

            switch(datatype)
            {
                case FLOAT:
                {
                    float norm, norm_A, eps, resid1, resid2;
                    eps = fla_lapack_slamch("P");

                    /* Test 1
                       compute norm(A - ((Q * lambda) * Q')) /
                                    (n * norm(A) * EPS)      */
                    norm_A = fla_lapack_slange("1", &n, &n, A, &lda, work);
                    sgemm_("N", "T", &n, &n, &n, &s_one, Q, &lda, Z, &lda, &s_n_one, A, &lda);
                    norm = fla_lapack_slange("1", &n, &n, A, &lda, work);
                    resid1 = norm / (eps * norm_A * (float)n);

                    /* Test 2
                       compute norm(I - Z'*Z) / (N * EPS)  */
                    resid2 = (float)check_orthogonality(datatype, Z, n, n, lda);

                    *residual = (double)fla_max(resid1, resid2);
                    break;
                }

                case DOUBLE:
                {
                    double norm, norm_A, eps, resid1, resid2;
                    eps = fla_lapack_dlamch("P");

                    /* Test 1
                       compute norm(A - (Q * lambda * Q')) /
                                     (n * norm(A) * EPS)       */
                    norm_A = fla_lapack_dlange("1", &n, &n, A, &lda, work);
                    dgemm_("N", "T", &n, &n, &n, &d_one, Q, &lda, Z, &lda, &d_n_one, A, &lda);
                    norm = fla_lapack_dlange("1", &n, &n, A, &lda, work);
                    resid1 = norm / (eps * norm_A * (double)n);

                    /* Test 2
                       compute norm(I - Z'*Z) / (N * EPS)*/
                    resid2 = check_orthogonality(datatype, Z, n, n, lda);

                    *residual = fla_max(resid1, resid2);
                    break;
                }

                case COMPLEX:
                {
                    float norm, norm_A, eps, resid1, resid2;
                    eps = fla_lapack_slamch("P");

                    /* Test 1
                       compute norm(A - (Q * lambda * Q')) /
                                       (n * norm(A) * EPS)   */
                    norm_A = fla_lapack_clange("1", &n, &n, A, &lda, work);
                    cgemm_("N", "C", &n, &n, &n, &c_one, Q, &lda, Z, &lda, &c_n_one, A, &lda);
                    norm = fla_lapack_clange("1", &n, &n, A, &lda, work);
                    resid1 = norm / (eps * norm_A * (float)n);

                    /* Test 2
                       compute norm(I - Z'*Z) / (N * EPS)*/
                    resid2 = (float)check_orthogonality(datatype, Z, n, n, lda);

                    *residual = (double)fla_max(resid1, resid2);
                    break;
                }

                case DOUBLE_COMPLEX:
                {
                    double norm, norm_A, eps, resid1, resid2;
                    eps = fla_lapack_dlamch("P");

                    /* Test 1
                       compute norm(A - (Q * lambda * Q')) /
                                    (n * norm(A) * EPS)      */
                    norm_A = fla_lapack_zlange("1", &n, &n, A, &lda, work);
                    zgemm_("N", "C", &n, &n, &n, &z_one, Q, &lda, Z, &lda, &z_n_one, A, &lda);
                    norm = fla_lapack_zlange("1", &n, &n, A, &lda, work);
                    resid1 = norm / (eps * norm_A * (double)n);

                    /* Test 2
                       compute norm(I - Z'*Z) / (N * EPS)*/
                    resid2 = check_orthogonality(datatype, Z, n, n, lda);

                    *residual = fla_max(resid1, resid2);
                    break;
                }
            }
            /* Test 3
               check if any of the eigen vectors failed to converge.
               Note: When eigen vectors fail to converge the corresponding
               element in ifail != 0 */
            if(ifail != NULL)
            {
                for(i = 0; i < n; i++)
                {
                    if(buff[i] != 0)
                        *residual = DBL_MAX;
                }
            }

            free_matrix(Q);
            free_matrix(Z);
        }
        /* Test 4: In case of specific input generation, compare input and
           output eigen values */
        if(L != NULL)
        {
            void *work = NULL;
            if(get_realtype(datatype) == FLOAT)
            {
                float norm, norm_L, eps, resid3;
                eps = fla_lapack_slamch("P");

                if((imatrix == 'O' || imatrix == 'U') && (scal != NULL))
                {
                    sscal_(&n, scal, L, &i_one);
                }

                norm_L = fla_lapack_slange("1", &n, &i_one, L, &i_one, work);
                saxpy_(&n, &s_n_one, lambda, &i_one, L, &i_one);
                norm = fla_lapack_slange("1", &n, &i_one, L, &i_one, work);
                resid3 = norm / (eps * norm_L * n);

                *residual = fla_max(*residual, (double)resid3);
            }
            else
            {
                double norm, norm_L, eps, resid3;
                eps = fla_lapack_dlamch("P");

                if((imatrix == 'O' || imatrix == 'U') && (scal != NULL))
                {
                    dscal_(&n, scal, L, &i_one);
                }
                norm_L = fla_lapack_dlange("1", &n, &i_one, L, &i_one, work);
                daxpy_(&n, &d_n_one, lambda, &i_one, L, &i_one);
                norm = fla_lapack_dlange("1", &n, &i_one, L, &i_one, work);
                resid3 = norm / (eps * norm_L * n);

                *residual = fla_max(*residual, resid3);
            }
        }
    }
}
