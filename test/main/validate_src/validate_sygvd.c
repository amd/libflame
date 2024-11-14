/******************************************************************************
 * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/* > \brief \b validate_sygvd.c                                              */
/* =========== DOCUMENTATION ===========                                     */
/* Definition:                                                               */
/* ===========                                                               */
/* SUBROUTINE  validate_sygvd(itype, &jobz, &range, n, A, A_test, lda,       */
/*                            B, B_test, ldb, Q, ldq, vl, vu, il, iu, L,     */
/*                            w, datatype, residual, imatrix, scal);         */
/* > \par Purpose:                                                           */
/* =============                                                             */
/* >                                                                         */
/* > \verbatim                                                               */
/* >                                                                         */
/* > Defines validate function of generalized symmetric-definite             */
/* > eigenproblem, of the form A*x=(lambda)*B*x,  A*Bx=(lambda)*x,           */
/* > or B*A*x=(lambda)*x                                                     */
/* >                                                                         */
/* > This API validates symmetric eigen APIs functionality by performing     */
/* > the folllowing tests:                                                   */
/* >                                                                         */
/* > TEST 1:                                                                 */
/* >    ITYPE 1: Verify A * X = (B * X * lambda)                             */
/* >    ITYPE 2: Verify A * B * X = X * lambda                               */
/* >    ITYPE 3: Verify B * A * X = (X * lambda)                             */
/* >             using norm calculations,                                    */
/* >                                                                         */
/* > TEST 2: Check if X * inv(X) = I                                         */
/* >                                                                         */
/* > TEST 3: Check B = LU, where L = U' are the calculated cholesky factors  */
/* >                                                                         */
/* > Test 4: check if any of the eigen vectors failed to converge.           */
/* >         Note: When eigen vectors fail to converge the corresponding     */
/* >               element in ifail != 0                                     */
/* >                                                                         */
/* > TEST 5: verify the input and output EVs are same or not                 */
/* >                                                                         */
/* > TEST I: When range = I, verify if EVs between the speficied indices     */
/* >         il, iu are generated                                            */
/* >                                                                         */
/* > \endverbatim                                                            */
/* ========================================================================= */

#include "test_common.h"

#define GET_TRANS_STR(datatype) (((datatype) == FLOAT || (datatype) == DOUBLE) ? "T" : "C")

#define invoke_lamch(type_prefix, arg) fla_lapack_##type_prefix##lamch(arg)
#define invoke_lange(type_prefix, ...) fla_lapack_##type_prefix##lange(__VA_ARGS__)
/* Invokes gemm for C = op(A) * op(B) - C */
#define invoke_gemm_diff(type_prefix, transA, transB, m, n, k, A, lda, B, ldb, C, ldc) \
    type_prefix##gemm_(transA, transB, m, n, k, &type_prefix##_one, A, lda, B, ldb,    \
                       &type_prefix##_n_one, C, ldc)
#define invoke_scal(type_prefix, ...) type_prefix##scal_(__VA_ARGS__)
#define invoke_axpy(type_prefix, ...) type_prefix##axpy_(__VA_ARGS__)

/* realtype: FLOAT/COMPLEX: float; DOUBLE, DOUBLE_COMPLEX: double
   realtype_prefix: FLOAT/COMPLEX: s; DOUBLE, DOUBLE_COMPLEX: d
   type_prefix: FLOAT: s; DOUBLE: d; COMPLEX: c; DOUBLE_COMPLEX: z
*/
#define test_1_body(realtype, realtype_prefix, type_prefix)                                 \
    {                                                                                       \
        realtype norm, norm_orig, resid;                                                    \
        /* Test 1 */                                                                        \
        /* Calculating LHS part of equation based on itype */                               \
        switch(itype)                                                                       \
        {                                                                                   \
            case 1:                                                                         \
                /* Z = A * X */                                                             \
                fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, A, &lda, X, &lda, Z, &lda); \
                break;                                                                      \
            case 2:                                                                         \
                /* Z = A * B * X */                                                         \
                fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, A, &lda, B, &ldb, P, &lda); \
                fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, P, &lda, X, &lda, Z, &lda); \
                break;                                                                      \
            case 3:                                                                         \
                /* Z = B * A * X */                                                         \
                fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, B, &ldb, A, &lda, P, &lda); \
                fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, P, &lda, X, &lda, Z, &lda); \
                break;                                                                      \
        }                                                                                   \
        norm_orig = invoke_lange(type_prefix, "1", &n, &n, Z, &lda, work);                  \
        if(norm_orig < ufmin)                                                               \
        {                                                                                   \
            norm_orig = ufmin;                                                              \
        }                                                                                   \
        /* F = X * lambda */                                                                \
        multiply_matrix_diag_vector(datatype, n, n, X, lda, lambda_out, 1);                 \
        switch(itype)                                                                       \
        {                                                                                   \
            case 1:                                                                         \
                /* P = B * F = B * X * lambda */                                            \
                fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, B, &ldb, X, &lda, P, &lda); \
                /* Z = Z - P = Z - (B * X * lambda ) */                                     \
                matrix_difference(datatype, n, n, Z, lda, P, lda);                          \
                break;                                                                      \
            case 2:                                                                         \
            case 3:                                                                         \
                /* Z = Z - (X * lambda) */                                                  \
                matrix_difference(datatype, n, n, Z, lda, X, lda);                          \
                break;                                                                      \
        }                                                                                   \
        norm = invoke_lange(type_prefix, "1", &n, &n, Z, &lda, work);                       \
        resid = norm / (eps * norm_orig * (realtype)n);                                     \
        *residual = (double)fla_max(*residual, resid);                                      \
    }

#define test_2_body(realtype, realtype_prefix, type_prefix)                                 \
    {                                                                                       \
        realtype norm, resid;                                                               \
        /* Test 2 */                                                                        \
        /* compute norm(X * inv(X) - I) / (N * EPS)  */                                     \
        /* Z = I */                                                                         \
        set_identity_matrix(datatype, n, n, Z, lda);                                        \
        /* Set X to the eigen values */                                                     \
        copy_matrix(datatype, "full", n, n, A_test, lda, X, lda);                           \
        /* Z = X * inv(X) - Z = X * inv(X) - I */                                           \
        invoke_gemm_diff(type_prefix, "N", "N", &n, &n, &n, X_inv, &lda, X, &lda, Z, &lda); \
        norm = invoke_lange(type_prefix, "1", &n, &n, Z, &lda, work);                       \
        resid = norm / (eps * (realtype)n);                                                 \
        *residual = (double)fla_max(*residual, resid);                                      \
    }

#define test_3_body(realtype, realtype_prefix, type_prefix)                             \
    {                                                                                   \
        realtype norm, resid;                                                           \
        /* Test 3 */                                                                    \
        /* Compute norm (LU - B) / (N * EPS * normB) */                                 \
        reset_matrix(datatype, n, n, Z, lda);                                           \
        copy_matrix(datatype, "full", n, n, B, ldb, Z, lda);                            \
        realtype normB = invoke_lange(type_prefix, "1", &n, &n, B, &lda, work);         \
        invoke_gemm_diff(type_prefix, "N", "N", &n, &n, &n, L, &lda, U, &lda, Z, &lda); \
        norm = invoke_lange(type_prefix, "1", &n, &n, Z, &lda, work);                   \
        resid = norm / (eps * normB * (realtype)n);                                     \
        *residual = (double)fla_max(*residual, resid);                                  \
    }

#define test_eigenvalues(realtype, realtype_prefix)                                         \
    {                                                                                       \
        realtype norm, norm_L, eps, resid3;                                                 \
        eps = invoke_lamch(realtype_prefix, "P");                                           \
        if(itype == 2 || itype == 3)                                                        \
        {                                                                                   \
            invoke_scal(realtype_prefix, &n, scal, lambda_orig, &i_one);                    \
        }                                                                                   \
        norm_L = invoke_lange(realtype_prefix, "1", &n, &i_one, lambda_orig, &i_one, work); \
        invoke_axpy(realtype_prefix, &n, &realtype_prefix##_n_one, lambda_out, &i_one,      \
                    lambda_orig, &i_one);                                                   \
        norm = invoke_lange(realtype_prefix, "1", &n, &i_one, lambda_orig, &i_one, work);   \
        resid3 = norm / (eps * norm_L * n);                                                 \
        *residual = fla_max(*residual, (double)resid3);                                     \
    }

#define invoke_tests(realtype, realtype_prefix, type_prefix) \
    {                                                        \
        realtype eps = invoke_lamch(realtype_prefix, "P");   \
        realtype ufmin = invoke_lamch(realtype_prefix, "U"); \
        test_1_body(realtype, realtype_prefix, type_prefix); \
        test_2_body(realtype, realtype_prefix, type_prefix); \
        test_3_body(realtype, realtype_prefix, type_prefix); \
    }

void validate_sygvd(integer itype, char *jobz, char *range, char *uplo, integer n, void *A,
                    void *A_test, integer lda, void *B, void *B_test, integer ldb, integer il,
                    integer iu, void *lambda_orig, void *lambda_out, void *ifail, integer datatype,
                    double *residual, char imatrix, void *scal)
{
    if(n == 0)
        return;
    *residual = 0;

    if(lambda_orig != NULL)
    {
        sort_realtype_vector(datatype, "A", n, lambda_orig, 1);
    }

    if((*range == 'I') || (*range == 'i'))
    {
        /* Test I range
           check if output EVs matches the input EVs in given index range */
        if((lambda_orig != NULL)
           && compare_realtype_vector(datatype, (iu - il + 1), lambda_orig, 1, il, lambda_out, 1))
        {
            *residual = DBL_MAX;
        }
    }
    else /* range A or V */
    {
        if(*jobz != 'N')
        {
            void *Z = NULL, *work = NULL, *X = NULL, *X_inv = NULL, *P = NULL, *L = NULL, *U = NULL;
            integer i;
            integer *buff = (integer *)ifail;

            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z, lda);
            reset_matrix(datatype, n, n, Z, lda);

            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &P, lda);
            reset_matrix(datatype, n, n, P, lda);

            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &X, lda);
            reset_matrix(datatype, n, n, X, lda);
            /* Copy the eigen vectors to X */
            copy_matrix(datatype, "full", n, n, A_test, lda, X, lda);

            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &X_inv, lda);
            reset_matrix(datatype, n, n, X_inv, lda);

            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &L, lda);
            reset_matrix(datatype, n, n, L, lda);

            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &U, lda);
            reset_matrix(datatype, n, n, U, lda);

            /* B = U'U = LL' = LU */
            if(*uplo == 'U')
            {
                /* B_test contains the upper triangular cholesky factor
                   Set L = U' */
                copy_matrix(datatype, "U", n, n, B_test, ldb, U, lda);
                set_identity_matrix(datatype, n, n, L, lda);
                /* L = U' * I */
                fla_invoke_trmm(datatype, "L", "U", GET_TRANS_STR(datatype), "N", &n, &n, U, &lda,
                                L, &lda);
            }
            else
            {
                /* B_test contains the lower triangular cholesky factor
                   Set U = L' */
                copy_matrix(datatype, "L", n, n, B_test, ldb, L, lda);
                set_identity_matrix(datatype, n, n, U, lda);
                /* U = L' * I */
                fla_invoke_trmm(datatype, "L", "L", GET_TRANS_STR(datatype), "N", &n, &n, L, &lda,
                                U, &lda);
            }

            /* Z = I */
            set_identity_matrix(datatype, n, n, Z, lda);
            /* P = X' */
            fla_invoke_gemm(datatype, GET_TRANS_STR(datatype), "N", &n, &n, &n, X, &lda, Z, &lda, P,
                            &lda);

            /* Getting X_inv */
            switch(itype)
            {
                case 1:
                case 2:
                    /* inv(X) = X' B = X' * L * U */
                    /* P = X' * L */
                    fla_invoke_trmm(datatype, "R", "L", "N", "N", &n, &n, L, &lda, P, &lda);
                    fla_invoke_trmm(datatype, "R", "U", "N", "N", &n, &n, U, &lda, P, &lda);
                    /* P = X' * L * U */
                    copy_matrix(datatype, "full", n, n, P, lda, X_inv, lda);
                    break;
                case 3:
                    /* inv(X) = X' inv(B) = X' inv(U) * inv(L) */
                    /* P = X' * inv(U) */
                    fla_invoke_trsm(datatype, "R", "U", "N", "N", &n, &n, U, &lda, P, &lda);
                    copy_matrix(datatype, "full", n, n, P, lda, X_inv, lda);
                    /* X_inv = P * inv(L) = X' * inv(U) * inv(L) */
                    fla_invoke_trsm(datatype, "R", "L", "N", "N", &n, &n, L, &lda, X_inv, &lda);
                    break;
            }

            switch(datatype)
            {
                case FLOAT:
                    invoke_tests(float, s, s);
                    break;
                case DOUBLE:
                    invoke_tests(double, d, d);
                    break;
                case COMPLEX:
                    invoke_tests(float, s, c);
                    break;
                case DOUBLE_COMPLEX:
                    invoke_tests(double, d, z);
                    break;
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
                        *residual = DBL_MAX;
                }
            }

            free_matrix(X);
            free_matrix(Z);
            free_matrix(X_inv);
            free_matrix(P);
            free_matrix(L);
            free_matrix(U);
        }
        /* Test 5: In case of specific input generation, compare input and
           output eigen values */
        if(lambda_orig != NULL)
        {
            void *work = NULL;
            if(get_realtype(datatype) == FLOAT)
            {
                test_eigenvalues(float, s);
            }
            else
            {
                test_eigenvalues(double, d);
            }
        }
    }
}
