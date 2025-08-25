/******************************************************************************
 * * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/
/*! @file fla_lapack_cholesky_small_kernels.h
 *  @brief Cholesky factorization kernels for small inputs.
 *  *  */

#include "FLAME.h"

#if FLA_ENABLE_AMD_OPT
#define FLA_POTRI_SMALL_1X1(uplo, n, buff_A, ldim_A, info) \
    /* Compute A⁻¹ = Ainvᵀ * Ainv = Ainv * Ainv */ \
    (buff_A)[0] = 1.0 / (((buff_A)[0]) * ((buff_A)[0])); \
    *(info) = 0;

#define FLA_POTRI_SMALL_2X2(datatype, uplo, n, buff_A, ldim_A, info) \
    char u = toupper(*(uplo)); \
    datatype invL[2][2] = {0}, *A = buff_A; \
    integer lda = *ldim_A; \
    \
    if (u == 'L') { \
        /* Invert lower triangular matrix L */ \
        invL[0][0] = 1.0 / A[0]; \
        invL[0][1] = 0.0; \
        invL[1][0] = -A[1] * invL[0][0] / A[(1 + lda)]; \
        invL[1][1] = 1.0 / A[1 + lda]; \
        \
        /* Compute A⁻¹ = invLᵀ * invL */ \
        A[0] = invL[0][0]*invL[0][0] + invL[1][0]*invL[1][0]; \
        A[1] = A[lda] = invL[0][0]*invL[0][1] + invL[1][0]*invL[1][1]; \
        A[1 + lda] = invL[0][1]*invL[0][1] + invL[1][1]*invL[1][1]; \
    } \
    else if (u == 'U') { \
        /* Invert upper triangular matrix U */ \
        invL[0][0] = 1.0 / A[0]; \
        invL[1][0] = 0.0; \
        invL[0][1] = -A[lda] * invL[0][0] / A[1 + lda]; \
        invL[1][1] = 1.0 / A[1 + lda]; \
        \
        /* Compute A⁻¹ = invU * invUᵀ */ \
        A[0] = invL[0][0]*invL[0][0] + invL[0][1]*invL[0][1]; \
        A[1] = A[lda] = invL[0][0]*invL[1][0] + invL[0][1]*invL[1][1]; \
        A[1 + lda] = invL[1][0]*invL[1][0] + invL[1][1]*invL[1][1]; \
    } \
    *(info) = 0;

#define FLA_POTRI_SMALL_3X3(datatype, uplo, n, buff_A, ldim_A, info) \
    char u = toupper(*(uplo)); \
    integer lda = *ldim_A; \
    datatype *A = buff_A; \
    \
    if (u == 'L') { \
        /* Step 1: Invert the lower triangular matrix L */ \
        datatype inv00 = 1.0 / A[0]; \
        datatype inv11 = 1.0 / A[1 + lda]; \
        datatype inv22 = 1.0 / A[2 + 2*lda]; \
        datatype inv10 = -A[1] * inv00 * inv11; \
        datatype inv21 = -A[2 + lda] * inv11 * inv22; \
        datatype inv20 = (-A[2] * inv00 - A[2 + lda] * inv10) * inv22; \
        \
        /* Step 2: Compute A⁻¹ = L_inv^T * L_inv and store result */ \
        A[0] = inv00*inv00 + inv10*inv10 + inv20*inv20; \
        A[lda] = A[1] = inv10*inv11 + inv20*inv21; \
        A[2*lda] = A[2] = inv20*inv22; \
        A[1 + lda] = inv11*inv11 + inv21*inv21; \
        A[1 + 2*lda] = A[2 + lda] = inv21*inv22; \
        A[2 + 2*lda] = inv22*inv22; \
    } \
    else if (u == 'U') { \
        /* Step 1: Invert the upper triangular matrix U */ \
        datatype inv00 = 1.0 / A[0]; \
        datatype inv11 = 1.0 / A[1 + lda]; \
        datatype inv22 = 1.0 / A[2 + 2*lda]; \
        datatype inv01 = -A[lda] * inv00 * inv11; \
        datatype inv12 = -A[1 + 2*lda] * inv11 * inv22; \
        datatype inv02 = (-A[2*lda] * inv22 - A[lda] * inv12) * inv00; \
        \
        /* Step 2: Compute A⁻¹ = U_inv * U_inv^T and store result */ \
        A[0] = inv00*inv00 + inv01*inv01 + inv02*inv02; \
        A[lda] = A[1] = inv01*inv11 + inv02*inv12; \
        A[2*lda] = A[2] = inv02*inv22; \
        A[1 + lda] = inv11*inv11 + inv12*inv12; \
        A[1 + 2*lda] = A[2 + lda] = inv12*inv22; \
        A[2 + 2*lda] = inv22*inv22; \
    } \
    *(info) = 0;

#define FLA_POTRI_SMALL_4X4(datatype, uplo, n, buff_A, ldim_A, info) \
    char u = toupper(*(uplo)); \
    integer lda = *ldim_A; \
    datatype *A = buff_A; \
    \
    if (u == 'L') { \
        /* Step 1: Invert the lower triangular matrix L */ \
        /* Diagonal elements */ \
        datatype inv00 = 1.0 / A[0]; \
        datatype inv11 = 1.0 / A[1 + lda]; \
        datatype inv22 = 1.0 / A[2 + 2*lda]; \
        datatype inv33 = 1.0 / A[3 + 3*lda]; \
        \
        /* Off-diagonal elements (row by row) */ \
        datatype inv10 = -A[1] * inv00 * inv11; \
        datatype inv20 = (-A[2] * inv00 - A[2 + lda] * inv10) * inv22; \
        datatype inv21 = -A[2 + lda] * inv11 * inv22; \
        datatype inv30 = (-A[3] * inv00 - A[3 + lda] * inv10 - A[3 + 2*lda] * inv20) * inv33; \
        datatype inv31 = (-A[3 + lda] * inv11 - A[3 + 2*lda] * inv21) * inv33; \
        datatype inv32 = -A[3 + 2*lda] * inv22 * inv33; \
        \
        /* Step 2: Compute A⁻¹ = L_inv^T * L_inv and store result */ \
        A[0] = inv00*inv00 + inv10*inv10 + inv20*inv20 + inv30*inv30; \
        A[lda] = A[1] = inv10*inv11 + inv20*inv21 + inv30*inv31; \
        A[2*lda] = A[2] = inv20*inv22 + inv30*inv32; \
        A[3*lda] = A[3] = inv30*inv33; \
        A[1 + lda] = inv11*inv11 + inv21*inv21 + inv31*inv31; \
        A[1 + 2*lda] = A[2 + lda] = inv21*inv22 + inv31*inv32; \
        A[1 + 3*lda] = A[3 + lda] = inv31*inv33; \
        A[2 + 2*lda] = inv22*inv22 + inv32*inv32; \
        A[2 + 3*lda] = A[3 + 2*lda] = inv32*inv33; \
        A[3 + 3*lda] = inv33*inv33; \
    } \
    else if (u == 'U') { \
        /* Step 1: Invert the upper triangular matrix U */ \
        /* Diagonal elements */ \
        datatype inv00 = 1.0 / A[0]; \
        datatype inv11 = 1.0 / A[1 + lda]; \
        datatype inv22 = 1.0 / A[2 + 2*lda]; \
        datatype inv33 = 1.0 / A[3 + 3*lda]; \
        \
        /* Off-diagonal elements (following exact 3x3 pattern) */ \
        datatype inv01 = -A[lda] * inv00 * inv11; \
        datatype inv12 = -A[1 + 2*lda] * inv11 * inv22; \
        datatype inv23 = -A[2 + 3*lda] * inv22 * inv33; \
        datatype inv02 = (-A[2*lda] * inv22 - A[lda] * inv12) * inv00; \
        datatype inv13 = (-A[1 + 3*lda] * inv33 - A[1 + 2*lda] * inv23) * inv11; \
        datatype inv03 = (-A[3*lda] * inv33 - A[lda] * inv13 - A[2*lda] * inv23) * inv00; \
        \
        /* Step 2: Compute A⁻¹ = U_inv * U_inv^T and store result */ \
        A[0] = inv00*inv00 + inv01*inv01 + inv02*inv02 + inv03*inv03; \
        A[lda] = A[1] = inv01*inv11 + inv02*inv12 + inv03*inv13; \
        A[2*lda] = A[2] = inv02*inv22 + inv03*inv23; \
        A[3*lda] = A[3] = inv03*inv33; \
        A[1 + lda] = inv11*inv11 + inv12*inv12 + inv13*inv13; \
        A[1 + 2*lda] = A[2 + lda] = inv12*inv22 + inv13*inv23; \
        A[1 + 3*lda] = A[3 + lda] = inv13*inv33; \
        A[2 + 2*lda] = inv22*inv22 + inv23*inv23; \
        A[2 + 3*lda] = A[3 + 2*lda] = inv23*inv33; \
        A[3 + 3*lda] = inv33*inv33; \
    } \
    *(info) = 0;

#endif // FLA_ENABLE_AMD_OPT