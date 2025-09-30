/*
    Copyright (C) 2025, Advanced Micro Devices, Inc.  All rights reserved.
*/
#include "test_bit_reproducibility.h"
#include "test_lapack.h"

/* Given a file path, create all constituent directories if missing */
void create_file_path_dirs(char *file_path)
{
    size_t len = strlen(file_path) + 1;
    char *tmp = malloc(len);
    if(file_path != NULL && tmp != NULL)
    {
        memcpy(tmp, file_path, len);
    }
    else
    {
        fprintf(stderr, "ERROR: Invalid file path\n");
        return;
    }
    char *p = tmp;
    while((p = strchr(p + 1, '/')) != NULL)
    {
        *p = 0;
        if(mkdir(tmp, S_IRWXU | S_IRWXG | S_IROTH) != 0 && errno != EEXIST)
        {
            fprintf(stderr, "Failed to create directory '%s': %s\n", tmp, strerror(errno));
        }
        *p = '/';
    }
    free(tmp);
}

/*
 * Allocates the required memory and generates a unique filename for the current test for
 * identification and running tests in batches. The function follows this flow:
 * 1. Iterate through all the params and Calculate the length of the filename
 * 2. Allocate memory for the filename
 * 3. Append all the strings and params to Generate the filename
 *
 * The File is created within folders for
 *     File type - "Ground-Truth" or "Input-Matrix", which refers to the type of file being created
 *     Test API - The name of the test API being run
 *     Seed - The seed used to generate the input data
 *     Data type - The data type of the input data
 *
 * Post this folder organization, each file is created with a unique name based on the parameters
 * passed to the API. This function supports variable number of input parameters for filename
 * creation. But for identification purposes, we need to send a string which informs the function
 * what is the type of parameter which is being passed. For the parameters, the following format
 * specifiers are supported: s - String d - Integer f - Float c - Char
 *
 *  */
void generate_filename(void **buffer, char *tst_api, int seed, char datatype, char imatrix_char,
                       const char *format, ...)
{

    /* Calculate the length of the filename */
    size_t buffer_len = 0;
    char seed_char[32];
    char datatype_char[2] = {(char)datatype, '\0'};

    snprintf(seed_char, sizeof(seed_char), "%d", seed);

    /* Along with the tst_api, seed and datatype, the 20 additional spaces are for the '/'
     * seperators and file type specifiers which will be added seperately */
    buffer_len = strlen(tst_api) + strlen(seed_char) + strlen(datatype_char) + 30;

    /* Iterating through all the params to collect required space for adding all params to filename
     */
    va_list args;
    va_start(args, format);
    for(const char *p = format; *p != '\0'; p++)
    {
        switch(*p)
        {
            case 's':
            {
                const char *str = va_arg(args, const char *);
                buffer_len = buffer_len + strlen(str);
                break;
            }
            case 'd':
            {
                int i = va_arg(args, int);
                char temp[32];
                snprintf(temp, sizeof(temp), "%d", i);
                buffer_len = buffer_len + strlen(temp);
                break;
            }
            case 'f':
            {
                double f = va_arg(args, double);
                char temp[64];
                snprintf(temp, sizeof(temp), "%.6g", f);
                buffer_len = buffer_len + strlen(temp);
                break;
            }
            case 'c':
            {
                int ch = va_arg(args, int);
                char temp[2] = {(char)ch, '\0'};
                buffer_len = buffer_len + strlen(temp);
                break;
            }
            default:
                // Ignore unknown format specifiers
                break;
        }
    }
    buffer_len = buffer_len + strlen(format) - 1; // For the '-' between the parameters
    buffer_len = buffer_len + 5; // For the '.bin' extension

    if(imatrix_char != '\0')
        buffer_len = buffer_len + 2; // For the '-' and the imatrix_char

    va_end(args);

    /* Allocate memory for the filename */
    *buffer = malloc(buffer_len * sizeof(char));

    va_start(args, format);

    /* The filetype, api name, seed and datatype are used for organizing the generated artifacts
     * into folders */
    strncpy((char *)(*buffer), tst_api, buffer_len);
    strncat((char *)(*buffer), "/", buffer_len - strlen((char *)(*buffer)) - 1);
    strncat((char *)(*buffer), seed_char, buffer_len - strlen((char *)(*buffer)) - 1);
    strncat((char *)(*buffer), "/", buffer_len - strlen((char *)(*buffer)) - 1);
    strncat((char *)(*buffer), datatype_char, buffer_len - strlen((char *)(*buffer)) - 1);
    strncat((char *)(*buffer), "/", buffer_len - strlen((char *)(*buffer)) - 1);
    /* Current filename "<filetype>/<tst_api>/<seed>/<datatype>/" */

    /* The params being passed to the API are used for creating the unique filename */
    int first = TRUE;
    for(const char *p = format; *p != '\0'; p++)
    {
        /* Add a '-' between each parameter */
        if(!first)
        {
            strncat((char *)(*buffer), "-", buffer_len - strlen((char *)(*buffer)) - 1);
        }

        switch(*p)
        {
            case 's':
            {
                /* Add the string param to the filename */
                const char *str = va_arg(args, const char *);
                strncat((char *)(*buffer), str, buffer_len - strlen((char *)(*buffer)) - 1);
                break;
            }
            case 'd':
            {
                /* Add the integer param to the filename */
                int i = va_arg(args, int);
                char temp[32];
                snprintf(temp, sizeof(temp), "%d", i);
                strncat((char *)(*buffer), temp, buffer_len - strlen((char *)(*buffer)) - 1);
                break;
            }
            case 'f':
            {
                /* Add the float param to the filename */
                double f = va_arg(args, double);
                char temp[64];
                snprintf(temp, sizeof(temp), "%.6g", f);
                strncat((char *)(*buffer), temp, buffer_len - strlen((char *)(*buffer)) - 1);
                break;
            }
            case 'c':
            {
                /* Add the char param to the filename */
                int ch = va_arg(args, int);
                char temp[2] = {(char)ch, '\0'};
                strncat((char *)(*buffer), temp, buffer_len - strlen((char *)(*buffer)) - 1);
                break;
            }
            default:
                // Ignore unknown format specifiers
                break;
        }
        first = FALSE;
    }

    va_end(args);

    /* Add the imatrix_char to the filename if present */
    if(imatrix_char != '\0')
    {
        char temp[2] = {imatrix_char, '\0'};
        strncat((char *)(*buffer), "-", buffer_len - strlen((char *)(*buffer)) - 1);
        strncat((char *)(*buffer), temp, buffer_len - strlen((char *)(*buffer)) - 1);
    }

    /* Add the file extension to the filename */
    strncat((char *)(*buffer), ".bin", buffer_len - strlen((char *)(*buffer)) - 1);
}

/* Helper function to append or overwrite filetype to filename
 * Note : This function assumes that the buffer has enough space to store filetype
 *        and that the filetypes being passed are of same length
 * "Input-Matrix" and "Ground-Truth" used by BRT are of same length, so this assumption is valid
 * */
void update_filetype(void **buffer, char *filetype, integer contains_filetype)
{
    if(contains_filetype)
    {
        /* Overwrite existing filetype with new filetype */
        memcpy((char *)(*buffer), filetype, strlen(filetype));
    }
    else
    {
        /* Move current filename  in memory to make space for file type */
        memmove((char *)(*buffer) + strlen(filetype), (char *)(*buffer),
                strlen((char *)(*buffer)) + 1);
        /* Append the filetype to start of filename*/
        memcpy((char *)(*buffer), filetype, strlen(filetype));
    }
}

/* Helper function to store and verify seed generated input between GT and V runs */
integer store_load_input_matrices(char BRT_char, char *filename, integer num_matrices, ...)
{
    if(same_char(BRT_char, 'G') || same_char(BRT_char, 'F'))
    {
        create_file_path_dirs(filename);

        /* Open the file to add ground truth data */
        FILE *gt_file = fopen(filename, "wb");
        if(gt_file == NULL)
        {
            printf("Error opening : %s\n", filename);
            perror("");
            return 0;
        }

        va_list args;
        va_start(args, num_matrices);
        for(int i = 0; i < num_matrices; i++)
        {
            integer datatype = va_arg(args, integer);
            integer m = va_arg(args, integer);
            integer n = va_arg(args, integer);
            void *A = va_arg(args, void *);
            integer lda = va_arg(args, integer);
            FLA_STORE_MATRIX(datatype, m, n, A, lda);
        }
        va_end(args);

        fclose(gt_file);
    }
    else if(same_char(BRT_char, 'V') || same_char(BRT_char, 'M'))
    {
        FILE *gt_file = fopen(filename, "rb");
        if(gt_file == NULL)
        {
            /* If no input data is found abort test as GT run has not been run*/
            printf("Please ensure initial reference/ground-truth file is available from earlier "
                   "runs\n");
            perror("");
            return 0;
        }
        va_list args;
        va_start(args, num_matrices);
        for(int i = 0; i < num_matrices; i++)
        {
            integer datatype = va_arg(args, integer);
            integer m = va_arg(args, integer);
            integer n = va_arg(args, integer);
            void *A = va_arg(args, void *);
            integer lda = va_arg(args, integer);
            FLA_LOAD_MATRIX(datatype, m, n, A, lda);
        }
        va_end(args);

        fclose(gt_file);
    }
    return 1;
}

/* This macro generates the multiplication factor required to interpret the memory blocks of
 * different datatypes to be of uint32_t memory block*/
#define FLA_GET_DATATYPE_FACTOR                          \
    int mul = 1;                                         \
    switch(datatype)                                     \
    {                                                    \
        case INTEGER:                                    \
        {                                                \
            mul = (sizeof(integer) / sizeof(uint32_t));  \
            break;                                       \
        }                                                \
        case FLOAT:                                      \
        {                                                \
            mul = (sizeof(float) / sizeof(uint32_t));    \
            break;                                       \
        }                                                \
        case DOUBLE:                                     \
        {                                                \
            mul = (sizeof(double) / sizeof(uint32_t));   \
            break;                                       \
        }                                                \
        case COMPLEX:                                    \
        {                                                \
            mul = (sizeof(scomplex) / sizeof(uint32_t));  \
            break;                                       \
        }                                                \
        case DOUBLE_COMPLEX:                             \
        {                                                \
            mul = (sizeof(dcomplex) / sizeof(uint32_t)); \
            break;                                       \
        }                                                \
    }

/* Compares the data present in memory for 2 vertices for bitwise equality */
int bitwise_compare_vector(integer datatype, integer m, void *A, void *B)
{
    FLA_GET_DATATYPE_FACTOR

    int temp;
    temp = memcmp(((uint32_t *)A), ((uint32_t *)B), m * mul * sizeof(uint32_t));

    return temp;
}

/* Compares the data present in memory for 2 matrices for bitwise equality */
int bitwise_compare_matrix(integer datatype, integer m, integer n, void *A, integer lda, void *B,
                           integer ldb)
{
    FLA_GET_DATATYPE_FACTOR

    /* Iterate and compare memory column by column */
    for(size_t i = 0; i < n; i++)
    {
        int temp;
        temp = memcmp(&((uint32_t *)A)[i * lda * mul], &((uint32_t *)B)[i * ldb * mul],
                      m * mul * sizeof(uint32_t));
        /* Break and return as soon as mismatching values are found */
        if(temp != 0)
            return temp;
    }
    return 0;
}

/* Compares the data present in memory for 2 matrices for bitwise equality, ignoring the first nb
 * elements */
int bitwise_compare_matrix_no_nb_diag(integer datatype, integer m, integer n, integer nb, void *A,
                                      integer lda, void *B, integer ldb)
{
    FLA_GET_DATATYPE_FACTOR

    /* Iterate and compare memory column by column */
    for(size_t i = 0; i < n; i++)
    {
        int temp;
        if(i < nb)
        {
            /* Compare elements before diagonal element */
            temp = memcmp(&((uint32_t *)A)[i * lda * mul], &((uint32_t *)B)[i * ldb * mul],
                          i * mul * sizeof(uint32_t));
            /* Break and return as soon as mismatching values are found */
            if(temp != 0)
                return temp;
            /* Compare elements after diagonal element */
            temp = memcmp(&((uint32_t *)A)[i * lda * mul + (i + 1) * mul],
                          &((uint32_t *)B)[i * ldb * mul + (i + 1) * mul],
                          (m - (i + 1)) * mul * sizeof(uint32_t));
        }
        else
        {
            /* Compare complete column */
            temp = memcmp(&((uint32_t *)A)[i * lda * mul], &((uint32_t *)B)[i * ldb * mul],
                          m * mul * sizeof(uint32_t));
        }
        /* Break and return as soon as mismatching values are found */
        if(temp != 0)
            return temp;
    }
    return 0;
}

/* Macro for binary division, this adds the MSB of a uint_32 to the right of the remainder and runs
 * XOR operation with the polynomial divisor followed by a left shift and does this for each bit in
 * the value to be divided */
#define BINARY_DIVISION_UINT32(a)                  \
    for(int z = 0; z < 32; z++)                    \
    {                                              \
        if((a & 0x80000000) == 0x80000000)         \
        {                                          \
            remainder = remainder | 1;             \
        }                                          \
        if((remainder & 0x80000000) == 0x80000000) \
        {                                          \
            remainder = remainder ^ polynomial;    \
        }                                          \
        a = a << 1;                                \
        remainder = remainder << 1;                \
    }

/*
 * CRC32_POLYNOMIAL is the polynomial used for CRC calculation in this implementation.
 * Value 0xAD0424F3 is chosen for bit reproducibility and has been chosen from the list of
 * polynomials which are studied by Philip Koopman from Carnegie Melon University.
 * https://users.ece.cmu.edu/~koopman/crc/
 */
#define CRC32_POLYNOMIAL 0xad0424f3

/* The Cyclic redundancy check, checksum generation follows a bit by bit modulo-2 binary division
 * algorithm. This function returns the checksum generated for a vector */
uint32_t generate_crc_vector(integer datatype, integer m, void *A)
{
    uint32_t remainder = 0;
    uint32_t polynomial = CRC32_POLYNOMIAL;

    FLA_GET_DATATYPE_FACTOR

    /* Load the initial 32 bits as first remainder and start binary division.
     * For each subsequent bit,
     * - Left shift the remainder to bring space for the new bit
     * - Add the new bit from next value to the remainder
     * - Left shift the next value and iteratively add the last bit for binary division */
    remainder = ((uint32_t *)A)[0];
    if((remainder & 0x80000000) == 0x80000000)
        remainder = remainder ^ polynomial;
    remainder = remainder << 1;

    for(size_t i = 1; i < (m * mul); i++)
    {
        uint32_t a = ((uint32_t *)A)[i];
        BINARY_DIVISION_UINT32(a);
    }
    return remainder;
}

/* The Cyclic redundancy check, checksum generation follows a bit by bit modulo-2 binary division
 * algorithm. This function returns the checksum generated for a matrix */
uint32_t generate_crc_matrix(integer datatype, integer m, integer n, void *A, integer lda)
{
    uint32_t remainder = 0;
    uint32_t polynomial = CRC32_POLYNOMIAL;

    FLA_GET_DATATYPE_FACTOR

    /* Load the initial 32 bits as first remainder and start binary division.
     * For each subsequent bit,
     * - Left shift the remainder to bring space for the new bit
     * - Add the new bit from next value to the remainder
     * - Left shift the next value and iteratively add the last bit for binary division */
    remainder = ((uint32_t *)A)[0];
    if((remainder & 0x80000000) == 0x80000000)
        remainder = remainder ^ polynomial;
    remainder = remainder << 1;

    int first = TRUE;

    for(size_t i = 0; i < n; i++)
    {
        for(size_t j = 0; j < (m * mul); j++)
        {
            /* Skip the first value as it has already been loaded into remainder */
            if(first)
            {
                first = FALSE;
                continue;
            }
            uint32_t a = ((uint32_t *)A)[(i * lda * mul) + j];
            BINARY_DIVISION_UINT32(a);
        }
    }
    return remainder;
}

/* The Cyclic redundancy check, checksum generation follows a bit by bit modulo-2 binary division
 * algorithm. This function returns the checksum generated for a matrix ignoring the first nb
 * diagonal elements */
uint32_t generate_crc_matrix_no_nb_diag(integer datatype, integer m, integer n, integer nb, void *A,
                                        integer lda)
{
    uint32_t remainder = 0;
    uint32_t polynomial = CRC32_POLYNOMIAL;

    FLA_GET_DATATYPE_FACTOR

    /* The first "mul" elements are skipped as we skip the 1st diagonal element of the original
     * matrix */
    /* Load the initial 32 bits as first remainder and start binary division.
     * For each subsequent bit,
     * - Left shift the remainder to bring space for the new bit
     * - Add the new bit from next value to the remainder
     * - Left shift the next value and iteratively add the last bit for binary division */
    remainder = ((uint32_t *)A)[mul];
    if((remainder & 0x80000000) == 0x80000000)
        remainder = remainder ^ polynomial;
    remainder = remainder << 1;

    int first = TRUE;

    for(size_t i = 0; i < n; i++)
    {
        for(size_t j = 0; j < (m * mul); j++)
        {
            /* Skipping the 1st diagonal element as well as the first element used for crc
             * generation */
            if(first)
            {
                if(j == mul)
                    first = FALSE;
                continue;
            }
            /* The first nb diagonals are skipped */
            if(i == (j / mul) && i < nb)
            {
                continue;
            }
            uint32_t a = ((uint32_t *)A)[(i * lda * mul) + j];
            BINARY_DIVISION_UINT32(a);
        }
    }
    return remainder;
}

/* Base function used to store the outputs of an API which are to be verified */
void store_outputs_base(void *filename, void *params, int nMatrix, int nVector, ...)
{
    /* Create and open a file for storing Ground truth*/
    FLA_OPEN_GT_FILE_STORE

    va_list args;
    va_start(args, nVector);
    for(int i = 0; i < nMatrix; i++)
    {
        integer datatype = va_arg(args, integer);
        integer m = va_arg(args, integer);
        integer n = va_arg(args, integer);
        void *A = va_arg(args, void *);
        integer lda = va_arg(args, integer);
        FLA_STORE_BRT_MATRIX(datatype, m, n, A, lda)
    }
    for(int i = 0; i < nVector; i++)
    {
        integer datatype = va_arg(args, integer);
        integer m = va_arg(args, integer);
        void *A = va_arg(args, void *);
        FLA_STORE_BRT_VECTOR(datatype, m, A)
    }
    va_end(args);

    FLA_CLOSE_GT_FILE_STORE
}

/* Base function used to load and verify bit exact match of the outputs of an API */
integer check_reproducibility_base(void *filename, void *params, int nMatrix, int nVector, ...)
{
    /* Open the file for reading Ground truth */
    FLA_OPEN_GT_FILE_READ

    va_list args;
    va_start(args, nVector);
    for(int i = 0; i < nMatrix; i++)
    {
        integer datatype = va_arg(args, integer);
        integer m = va_arg(args, integer);
        integer n = va_arg(args, integer);
        void *A = va_arg(args, void *);
        integer lda = va_arg(args, integer);
        FLA_VERIFY_BRT_MATRIX(datatype, m, n, A, lda)
    }
    for(int i = 0; i < nVector; i++)
    {
        integer datatype = va_arg(args, integer);
        integer m = va_arg(args, integer);
        void *A = va_arg(args, void *);
        FLA_VERIFY_BRT_VECTOR(datatype, m, A)
    }
    va_end(args);

    fclose(gt_file);
    return 1;
}