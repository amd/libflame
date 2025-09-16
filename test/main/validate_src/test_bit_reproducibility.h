/******************************************************************************
 * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file test_bit_reproducibility.h
 *  @brief Defines function declarations for bit reproducibility test to use in APIs of test
 * suite.
 *  */

#ifndef TEST_BIT_REPRODUCIBILITY_H
#define TEST_BIT_REPRODUCIBILITY_H

#include "test_common.h"

/*
 * Windows-specific code block:
 * - Includes <Windows.h> and <direct.h> for Windows API and directory functions.
 * - Redefines mkdir to _mkdir for compatibility with POSIX code.
 * This is required because Windows does not provide the POSIX mkdir with two arguments.
 * If building on Windows (_WIN32 is defined), use these headers and macro.
 */
#ifdef _WIN32
#include <Windows.h>
#include <direct.h>
#define mkdir(dir, mode) _mkdir(dir)
#endif

/* Stores the complete binary data of a vector into an open "gt_file" */
#define FLA_STORE_VECTOR(datatype, m, A)                               \
    switch(datatype)                                                   \
    {                                                                  \
        case INTEGER:                                                  \
        {                                                              \
            fwrite(&((integer *)A)[0], sizeof(integer), m, gt_file);   \
            break;                                                     \
        }                                                              \
        case FLOAT:                                                    \
        {                                                              \
            fwrite(&((float *)A)[0], sizeof(float), m, gt_file);       \
            break;                                                     \
        }                                                              \
        case DOUBLE:                                                   \
        {                                                              \
            fwrite(&((double *)A)[0], sizeof(double), m, gt_file);     \
            break;                                                     \
        }                                                              \
        case COMPLEX:                                                  \
        {                                                              \
            fwrite(&((scomplex *)A)[0], sizeof(scomplex), m, gt_file); \
            break;                                                     \
        }                                                              \
        case DOUBLE_COMPLEX:                                           \
        {                                                              \
            fwrite(&((dcomplex *)A)[0], sizeof(dcomplex), m, gt_file); \
            break;                                                     \
        }                                                              \
    }

/* Loads the complete binary data of a vector from an open "gt_file" */
#define FLA_LOAD_VECTOR(datatype, m, A)                                          \
    switch(datatype)                                                             \
    {                                                                            \
        case INTEGER:                                                            \
        {                                                                        \
            if(!(fread(&((integer *)A)[0], sizeof(integer), m, gt_file) == m))   \
            {                                                                    \
                printf("Error: Incomplete data stored in GT_file\n");            \
                break;                                                           \
            }                                                                    \
            break;                                                               \
        }                                                                        \
        case FLOAT:                                                              \
        {                                                                        \
            if(!(fread(&((float *)A)[0], sizeof(float), m, gt_file) == m))       \
            {                                                                    \
                printf("Error: Incomplete data stored in GT_file\n");            \
                break;                                                           \
            }                                                                    \
            break;                                                               \
        }                                                                        \
        case DOUBLE:                                                             \
        {                                                                        \
            if(!(fread(&((double *)A)[0], sizeof(double), m, gt_file) == m))     \
            {                                                                    \
                printf("Error: Incomplete data stored in GT_file\n");            \
                break;                                                           \
            }                                                                    \
            break;                                                               \
        }                                                                        \
        case COMPLEX:                                                            \
        {                                                                        \
            if(!(fread(&((scomplex *)A)[0], sizeof(scomplex), m, gt_file) == m)) \
            {                                                                    \
                printf("Error: Incomplete data stored in GT_file\n");            \
                break;                                                           \
            }                                                                    \
            break;                                                               \
        }                                                                        \
        case DOUBLE_COMPLEX:                                                     \
        {                                                                        \
            if(!(fread(&((dcomplex *)A)[0], sizeof(dcomplex), m, gt_file) == m)) \
            {                                                                    \
                printf("Error: Incomplete data stored in GT_file\n");            \
                break;                                                           \
            }                                                                    \
            break;                                                               \
        }                                                                        \
    }

/* Stores the complete binary data of a matrix into an open "gt_file" */
#define FLA_STORE_MATRIX(datatype, m, n, A, lda)                                 \
    switch(datatype)                                                             \
    {                                                                            \
        case INTEGER:                                                            \
        {                                                                        \
            for(size_t i = 0; i < n; i++)                                        \
            {                                                                    \
                fwrite(&((integer *)A)[i * lda], sizeof(integer), m, gt_file);   \
            }                                                                    \
            break;                                                               \
        }                                                                        \
        case FLOAT:                                                              \
        {                                                                        \
            for(size_t i = 0; i < n; i++)                                        \
            {                                                                    \
                fwrite(&((float *)A)[i * lda], sizeof(float), m, gt_file);       \
            }                                                                    \
            break;                                                               \
        }                                                                        \
        case DOUBLE:                                                             \
        {                                                                        \
            for(size_t i = 0; i < n; i++)                                        \
            {                                                                    \
                fwrite(&((double *)A)[i * lda], sizeof(double), m, gt_file);     \
            }                                                                    \
            break;                                                               \
        }                                                                        \
        case COMPLEX:                                                            \
        {                                                                        \
            for(size_t i = 0; i < n; i++)                                        \
            {                                                                    \
                fwrite(&((scomplex *)A)[i * lda], sizeof(scomplex), m, gt_file); \
            }                                                                    \
            break;                                                               \
        }                                                                        \
        case DOUBLE_COMPLEX:                                                     \
        {                                                                        \
            for(size_t i = 0; i < n; i++)                                        \
            {                                                                    \
                fwrite(&((dcomplex *)A)[i * lda], sizeof(dcomplex), m, gt_file); \
            }                                                                    \
            break;                                                               \
        }                                                                        \
    }

/* Loads the complete binary data of a matrix from an open "gt_file" */
#define FLA_LOAD_MATRIX(datatype, m, n, A, lda)                                            \
    switch(datatype)                                                                       \
    {                                                                                      \
        case INTEGER:                                                                      \
        {                                                                                  \
            for(size_t i = 0; i < n; i++)                                                  \
            {                                                                              \
                if(!(fread(&((integer *)A)[i * lda], sizeof(integer), m, gt_file) == m))   \
                {                                                                          \
                    printf("Error: Incomplete data stored in GT_file\n");                  \
                    break;                                                                 \
                }                                                                          \
            }                                                                              \
            break;                                                                         \
        }                                                                                  \
        case FLOAT:                                                                        \
        {                                                                                  \
            for(size_t i = 0; i < n; i++)                                                  \
            {                                                                              \
                if(!(fread(&((float *)A)[i * lda], sizeof(float), m, gt_file) == m))       \
                {                                                                          \
                    printf("Error: Incomplete data stored in GT_file\n");                  \
                    break;                                                                 \
                }                                                                          \
            }                                                                              \
            break;                                                                         \
        }                                                                                  \
        case DOUBLE:                                                                       \
        {                                                                                  \
            for(size_t i = 0; i < n; i++)                                                  \
            {                                                                              \
                if(!(fread(&((double *)A)[i * lda], sizeof(double), m, gt_file) == m))     \
                {                                                                          \
                    printf("Error: Incomplete data stored in GT_file\n");                  \
                    break;                                                                 \
                }                                                                          \
            }                                                                              \
            break;                                                                         \
        }                                                                                  \
        case COMPLEX:                                                                      \
        {                                                                                  \
            for(size_t i = 0; i < n; i++)                                                  \
            {                                                                              \
                if(!(fread(&((scomplex *)A)[i * lda], sizeof(scomplex), m, gt_file) == m)) \
                {                                                                          \
                    printf("Error: Incomplete data stored in GT_file\n");                  \
                    break;                                                                 \
                }                                                                          \
            }                                                                              \
            break;                                                                         \
        }                                                                                  \
        case DOUBLE_COMPLEX:                                                               \
        {                                                                                  \
            for(size_t i = 0; i < n; i++)                                                  \
            {                                                                              \
                if(!(fread(&((dcomplex *)A)[i * lda], sizeof(dcomplex), m, gt_file) == m)) \
                {                                                                          \
                    printf("Error: Incomplete data stored in GT_file\n");                  \
                    break;                                                                 \
                }                                                                          \
            }                                                                              \
            break;                                                                         \
        }                                                                                  \
    }

/* Helper macro used in store_API_outputs function to store the CRC and binary data of a matrix */
#define FLA_STORE_BRT_MATRIX(datatype, m, n, A, lda)                   \
    if(same_char(((test_params_t *)params)->BRT_char, 'G'))            \
    {                                                                  \
        uint32_t A_GT;                                                 \
        A_GT = generate_crc_matrix(datatype, m, n, A, lda);            \
        fwrite(&A_GT, sizeof(uint32_t), 1, gt_file);                   \
    }                                                                  \
    else if(same_char(((test_params_t *)params)->BRT_char, 'F'))       \
    {                                                                  \
        FLA_STORE_MATRIX(datatype, m, n, A, lda);                      \
    }                                                                  \
    else if(same_char(((test_params_t *)params)->BRT_char, 'L'))       \
    {                                                                  \
        uint32_t A_GT;                                                 \
        A_GT = generate_crc_vector(datatype, m, A);                    \
        printf("Output: %dx%d matrix\tCRC:%" PRIu32 "\n", m, n, A_GT); \
    }

/* Helper macro used in store_API_outputs function to store the CRC and binary data of a matrix
 * without the first nb diagonal elements */
#define FLA_STORE_BRT_MATRIX_NB_DIAG(datatype, m, n, nb, A, lda)           \
    if(same_char(((test_params_t *)params)->BRT_char, 'G'))                \
    {                                                                      \
        uint32_t A_GT;                                                     \
        A_GT = generate_crc_matrix_no_nb_diag(datatype, m, n, nb, A, lda); \
        fwrite(&A_GT, sizeof(uint32_t), 1, gt_file);                       \
    }                                                                      \
    else if(same_char(((test_params_t *)params)->BRT_char, 'F'))           \
    {                                                                      \
        FLA_STORE_MATRIX(datatype, m, n, A, lda);                          \
    }                                                                      \
    else if(same_char(((test_params_t *)params)->BRT_char, 'L'))           \
    {                                                                      \
        uint32_t A_GT;                                                     \
        A_GT = generate_crc_matrix_no_nb_diag(datatype, m, n, nb, A, lda); \
        printf("Output: %dx%d matrix\tCRC:%" PRIu32 "\n", m, n, A_GT);     \
    }

/* Helper macro used in store_API_outputs function to store the CRC and binary data of a vector */
#define FLA_STORE_BRT_VECTOR(datatype, m, A)                     \
    if(same_char(((test_params_t *)params)->BRT_char, 'G'))      \
    {                                                            \
        uint32_t A_GT;                                           \
        A_GT = generate_crc_vector(datatype, m, A);              \
        fwrite(&A_GT, sizeof(uint32_t), 1, gt_file);             \
    }                                                            \
    else if(same_char(((test_params_t *)params)->BRT_char, 'F')) \
    {                                                            \
        FLA_STORE_VECTOR(datatype, m, A);                        \
    }                                                            \
    else if(same_char(((test_params_t *)params)->BRT_char, 'L')) \
    {                                                            \
        uint32_t A_GT;                                           \
        A_GT = generate_crc_vector(datatype, m, A);              \
        printf("Output: %d vector\tCRC:%" PRIu32 "\n", m, A_GT); \
    }

/* Helper macro used in check_bit_reproducibility_API function to verify the CRC and binary data of
 * a matrix */
#define FLA_VERIFY_BRT_MATRIX(datatype, m, n, A, lda)                      \
    if(same_char(((test_params_t *)params)->BRT_char, 'V'))                \
    {                                                                      \
        uint32_t A_GT;                                                     \
        if(!(fread(&A_GT, sizeof(uint32_t), 1, gt_file) == 1))             \
        {                                                                  \
            printf("Error: Incomplete data stored in GT_file\n");          \
            fclose(gt_file);                                               \
            return 0;                                                      \
        }                                                                  \
        if(A_GT != generate_crc_matrix(datatype, m, n, A, lda))            \
        {                                                                  \
            fclose(gt_file);                                               \
            return 0;                                                      \
        }                                                                  \
    }                                                                      \
    else if(same_char(((test_params_t *)params)->BRT_char, 'M'))           \
    {                                                                      \
        void *A_GT = NULL;                                                 \
        create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_GT, lda);       \
        FLA_LOAD_MATRIX(datatype, m, n, A_GT, lda)                         \
        if(bitwise_compare_matrix(datatype, m, n, A, lda, A_GT, lda) != 0) \
        {                                                                  \
            free_matrix(A_GT);                                             \
            fclose(gt_file);                                               \
            return 0;                                                      \
        }                                                                  \
        free_matrix(A_GT);                                                 \
    }

/* Helper macro used in check_bit_reproducibility_API function to verify the CRC and binary data of
 * a vector */
#define FLA_VERIFY_BRT_VECTOR(datatype, m, A)                     \
    if(same_char(((test_params_t *)params)->BRT_char, 'V'))       \
    {                                                             \
        uint32_t A_GT;                                            \
        if(!(fread(&A_GT, sizeof(uint32_t), 1, gt_file) == 1))    \
        {                                                         \
            printf("Error: Incomplete data stored in GT_file\n"); \
            fclose(gt_file);                                      \
            return 0;                                             \
        }                                                         \
        if(A_GT != generate_crc_vector(datatype, m, A))           \
        {                                                         \
            fclose(gt_file);                                      \
            return 0;                                             \
        }                                                         \
    }                                                             \
    else if(same_char(((test_params_t *)params)->BRT_char, 'M'))  \
    {                                                             \
        void *A_GT = NULL;                                        \
        create_vector(datatype, &A_GT, m);                        \
        FLA_LOAD_VECTOR(datatype, m, A_GT)                        \
        if(bitwise_compare_vector(datatype, m, A, A_GT) != 0)     \
        {                                                         \
            free_vector(A_GT);                                    \
            fclose(gt_file);                                      \
            return 0;                                             \
        }                                                         \
        free_vector(A_GT);                                        \
    }

/* Helper macro used in check_bit_reproducibility_API function to verify the CRC and binary data of
 * a matrix without the first nb diagonal elements */
#define FLA_VERIFY_BRT_MATRIX_NB_DIAG(datatype, m, n, nb, A, lda)                         \
    if(same_char(((test_params_t *)params)->BRT_char, 'V'))                               \
    {                                                                                     \
        uint32_t A_GT;                                                                    \
        if(!(fread(&A_GT, sizeof(uint32_t), 1, gt_file) == 1))                            \
        {                                                                                 \
            printf("Error: Incomplete data stored in GT_file\n");                         \
            fclose(gt_file);                                                              \
            return 0;                                                                     \
        }                                                                                 \
        if(A_GT != generate_crc_matrix_no_nb_diag(datatype, m, n, nb, A, lda))            \
        {                                                                                 \
            fclose(gt_file);                                                              \
            return 0;                                                                     \
        }                                                                                 \
    }                                                                                     \
    else if(same_char(((test_params_t *)params)->BRT_char, 'M'))                          \
    {                                                                                     \
        void *A_GT = NULL;                                                                \
        create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_GT, lda);                      \
        FLA_LOAD_MATRIX(datatype, m, n, A_GT, lda)                                        \
        if(bitwise_compare_matrix_no_nb_diag(datatype, m, n, nb, A, lda, A_GT, lda) != 0) \
        {                                                                                 \
            free_matrix(A_GT);                                                            \
            fclose(gt_file);                                                              \
            return 0;                                                                     \
        }                                                                                 \
        free_matrix(A_GT);                                                                \
    }

/* Macro for the BRT path in the post processing of API outputs in the test driver */
#define IF_FLA_BRT_VALIDATION(m, n, storeFunction, validateFunction, checkBRTFunction) \
    if(FLA_BIT_REPRODUCIBILITY_TEST)                                                   \
    {                                                                                  \
        if(same_char(params->BRT_char, 'G') || same_char(params->BRT_char, 'F')        \
           || same_char(params->BRT_char, 'L'))                                        \
        {                                                                              \
            storeFunction;                                                             \
            validateFunction;                                                          \
        }                                                                              \
        else if(same_char(params->BRT_char, 'V') || same_char(params->BRT_char, 'M'))  \
        {                                                                              \
            if(!checkBRTFunction)                                                      \
            {                                                                          \
                residual = DBL_MAX;                                                    \
            }                                                                          \
            else                                                                       \
            {                                                                          \
                residual = 0;                                                          \
            }                                                                          \
            FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);                         \
        }                                                                              \
    }

/* Macro for the BRT path in the pre processing of API inputs in the test driver
 * This macro generates the unique filename and stores the API input in GT case and loads API input
 * in V case */
#define FLA_BRT_PROCESS_SINGLE_INPUT(datatype, m, n, A, lda, format, ...)                         \
    if(FLA_BIT_REPRODUCIBILITY_TEST)                                                              \
    {                                                                                             \
        generate_filename(&filename, tst_api, params->seed, get_datatype_char(datatype),          \
                          params->imatrix_char, format, __VA_ARGS__);                             \
        if(!same_char(((test_params_t *)params)->BRT_char, 'L'))                                  \
        {                                                                                         \
            update_filetype(&filename, "BRT/Input-Matrix/", FALSE);                               \
            if(!store_load_input_matrices(params->BRT_char, filename, 1, datatype, m, n, A, lda)) \
            {                                                                                     \
                goto free_buffers;                                                                \
            }                                                                                     \
        }                                                                                         \
        else                                                                                      \
        {                                                                                         \
            ((char *)filename)[strlen(filename) - 4] = '\0';                                      \
        }                                                                                         \
    }

/* Macro for the BRT path in the pre processing of API inputs in the test driver
 * This macro generates the unique filename and stores the API input in GT case and loads API input
 * in V case */
#define FLA_BRT_PROCESS_TWO_INPUT(datatype, m, n, A, lda, datatypeB, mB, nB, B, ldb, format, ...) \
    if(FLA_BIT_REPRODUCIBILITY_TEST)                                                              \
    {                                                                                             \
        generate_filename(&filename, tst_api, params->seed, get_datatype_char(datatype),          \
                          params->imatrix_char, format, __VA_ARGS__);                             \
        if(!same_char(((test_params_t *)params)->BRT_char, 'L'))                                  \
        {                                                                                         \
            update_filetype(&filename, "BRT/Input-Matrix/", FALSE);                               \
            if(!store_load_input_matrices(params->BRT_char, filename, 2, datatype, m, n, A, lda,  \
                                          datatypeB, mB, nB, B, ldb))                             \
            {                                                                                     \
                goto free_buffers;                                                                \
            }                                                                                     \
        }                                                                                         \
        else                                                                                      \
        {                                                                                         \
            ((char *)filename)[strlen(filename) - 4] = '\0';                                      \
        }                                                                                         \
    }

/* Macro for the BRT path in the pre processing of API inputs in the test driver
 * This macro generates the unique filename and stores the API input in GT case and loads API input
 * in V case */
#define FLA_BRT_PROCESS_THREE_INPUT(datatype, m, n, A, lda, datatypeB, mB, nB, B, ldb, datatypeC, \
                                    mC, nC, C, ldc, format, ...)                                  \
    if(FLA_BIT_REPRODUCIBILITY_TEST)                                                              \
    {                                                                                             \
        generate_filename(&filename, tst_api, params->seed, get_datatype_char(datatype),          \
                          params->imatrix_char, format, __VA_ARGS__);                             \
        if(!same_char(((test_params_t *)params)->BRT_char, 'L'))                                  \
        {                                                                                         \
            update_filetype(&filename, "BRT/Input-Matrix/", FALSE);                               \
            if(!store_load_input_matrices(params->BRT_char, filename, 3, datatype, m, n, A, lda,  \
                                          datatypeB, mB, nB, B, ldb, datatypeC, mC, nC, C, ldc))  \
            {                                                                                     \
                goto free_buffers;                                                                \
            }                                                                                     \
        }                                                                                         \
        else                                                                                      \
        {                                                                                         \
            ((char *)filename)[strlen(filename) - 4] = '\0';                                      \
        }                                                                                         \
    }

/* Macro for the BRT path in the pre processing of API inputs in the test driver
 * This macro generates the unique filename and stores the API input in GT case and loads API input
 * in V case */
#define FLA_BRT_PROCESS_FOUR_INPUT(datatype, m, n, A, lda, datatypeB, mB, nB, B, ldb, datatypeC, \
                                   mC, nC, C, ldc, datatypeD, mD, nD, D, ldd, format, ...)       \
    if(FLA_BIT_REPRODUCIBILITY_TEST)                                                             \
    {                                                                                            \
        generate_filename(&filename, tst_api, params->seed, get_datatype_char(datatype),         \
                          params->imatrix_char, format, __VA_ARGS__);                            \
        if(!same_char(((test_params_t *)params)->BRT_char, 'L'))                                 \
        {                                                                                        \
            update_filetype(&filename, "BRT/Input-Matrix/", FALSE);                              \
            if(!store_load_input_matrices(params->BRT_char, filename, 4, datatype, m, n, A, lda, \
                                          datatypeB, mB, nB, B, ldb, datatypeC, mC, nC, C, ldc,  \
                                          datatypeD, mD, nD, D, ldd))                            \
            {                                                                                    \
                goto free_buffers;                                                               \
            }                                                                                    \
        }                                                                                        \
        else                                                                                     \
        {                                                                                        \
            ((char *)filename)[strlen(filename) - 5] = '\0';                                     \
        }                                                                                        \
    }

/* Free filename for BRT path*/
#define FLA_FREE_FILENAME(filename) \
    if(filename != NULL)            \
    {                               \
        free(filename);             \
        filename = NULL;            \
    }

/* Macro for opening file for writing ground truth */
#define FLA_OPEN_GT_FILE_STORE                                 \
    FILE *gt_file = NULL;                                      \
    if(!same_char(((test_params_t *)params)->BRT_char, 'L'))   \
    {                                                          \
        update_filetype(&filename, "BRT/Ground-Truth/", TRUE); \
        create_file_path_dirs((char *)filename);               \
        gt_file = fopen((char *)filename, "wb");               \
        if(gt_file == NULL)                                    \
        {                                                      \
            printf("Error opening : %s\n", (char *)filename);  \
            perror("");                                        \
            return;                                            \
        }                                                      \
    }                                                          \
    else                                                       \
    {                                                          \
        printf("Test ID: %s\n", (char *)filename);             \
    }

/* Macro for closing file after writing ground truth */
#define FLA_CLOSE_GT_FILE_STORE                              \
    if(!same_char(((test_params_t *)params)->BRT_char, 'L')) \
    {                                                        \
        fclose(gt_file);                                     \
    }

/* Macro for opening file for reading ground truth */
#define FLA_OPEN_GT_FILE_READ                              \
    update_filetype(&filename, "BRT/Ground-Truth/", TRUE); \
    FILE *gt_file = fopen(filename, "rb");                 \
    if(gt_file == NULL)                                    \
    {                                                      \
        printf("Error opening : %s\n", (char *)filename);  \
        perror("");                                        \
        return 0;                                          \
    }

/* Given a file path, create all constituent directories if missing */
void create_file_path_dirs(char *file_path);

/* Generates a unique filename for the current test for identification and running tests in batches.
 */
void generate_filename(void **buffer, char *tst_api, int seed, char datatype, char imatrix_char,
                       const char *format, ...);

/* Overwrites or appends filetype to the filename */
void update_filetype(void **buffer, char *filetype, integer contains_filetype);

/* Helper function to store and verify seed generated input between GT and V runs */
integer store_load_input_matrices(char BRT_char, char *filename, integer num_matrices, ...);

/* Compares the data present in memory for 2 vectors for bitwise equality */
int bitwise_compare_vector(integer datatype, integer m, void *A, void *B);

/* Compares the data present in memory for 2 matrices for bitwise equality */
int bitwise_compare_matrix(integer datatype, integer m, integer n, void *A, integer lda, void *B,
                           integer ldb);

/* Compares the data present in memory for 2 matrices for bitwise equality, ignoring the first nb
 * elements */
int bitwise_compare_matrix_no_nb_diag(integer datatype, integer m, integer n, integer nb, void *A,
                                      integer lda, void *B, integer ldb);

/* Generate CRC remainder for a vector */
uint32_t generate_crc_vector(integer datatype, integer m, void *A);

/* Generate CRC remainder for a matrix */
uint32_t generate_crc_matrix(integer datatype, integer m, integer n, void *A, integer lda);

/* Generate CRC remainder for a matrix excluding first nb values of the diagonal */
uint32_t generate_crc_matrix_no_nb_diag(integer datatype, integer m, integer n, integer nb, void *A,
                                        integer lda);

/* Base function used to store the outputs of an API which are to be verified
 *
 * nMatrix and nVector arguments specify the number of matrices and vectors which are passed to the
 * function. The following arguments are required in the same order for each matrix
 *    - datatype, rows, columns, buffer, leading dimension
 * Followed by the following arguments for each vector
 *    - datatype, number of elements, buffer
 *
 * All matrices must be defined before moving to vectors
 *
 * Note : The order of buffers must be the same for both store_outputs_base function and
 * check_reproducibility_base function
 * */
void store_outputs_base(void *filename, void *params, int nMatrix, int nVector, ...);

/* Base function used to load and verify bit exact match of the outputs of an API
 *
 * nMatrix and nVector arguments specify the number of matrices and vectors which are passed to the
 * function. The following arguments are required in the same order for each matrix
 *    - datatype, rows, columns, buffer, leading dimension
 * Followed by the following arguments for each vector
 *    - datatype, number of elements, buffer
 *
 * All matrices must be defined before moving to vectors
 *
 * Note : The order of buffers must be the same for both store_outputs_base function and
 * check_reproducibility_base function
 *  */
integer check_reproducibility_base(void *filename, void *params, int nMatrix, int nVector, ...);
#endif