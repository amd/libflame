/*
    Copyright (C) 2025, Advanced Micro Devices, Inc.  All rights reserved.
*/
#include "test_bit_reproducibility.h"
#include "test_lapack.h"

/* Given a file path, create all constituent directories if missing */
void create_file_path_dirs(char *file_path)
{
    char *tmp = strdup(file_path);
    if(tmp == NULL)
    {
        // Handle allocation failure
        fprintf(stderr, "Error: strdup failed in create_file_path_dirs\n");
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
    strncpy((char *)(*buffer), "/", buffer_len);
    strncat((char *)(*buffer), tst_api, buffer_len - strlen((char *)(*buffer)) - 1);
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

/* Counts the bits which are mismatched in 2 integer values from integer buffers
 *    - Load data into uint32_t to enable bitwise operations
 *    - Run XOR operation on the 2 uint32_t values for setting mismatch bits
 *    - Count the number of set bits in the XOR result to get the number of mismatch
 *  */
#define bitwise_compare_int(indexA, indexB)                   \
    {                                                         \
        uint32_t a, b;                                        \
        memcpy(&a, &((integer *)A)[indexA], sizeof(integer)); \
        memcpy(&b, &((integer *)B)[indexB], sizeof(integer)); \
        mismatchBits = mismatchBits + count_set_bits(a ^ b);  \
    }

/* Counts the bits which are mismatched in 2 float values from float buffers
 *    - Load data into uint32_t to enable bitwise operations
 *    - Run XOR operation on the 2 uint32_t values for setting mismatch bits
 *    - Count the number of set bits in the XOR result to get the number of mismatch
 *  */
#define bitwise_compare_float(indexA, indexB)                \
    {                                                        \
        uint32_t a, b;                                       \
        memcpy(&a, &((float *)A)[indexA], sizeof(float));    \
        memcpy(&b, &((float *)B)[indexB], sizeof(float));    \
        mismatchBits = mismatchBits + count_set_bits(a ^ b); \
    }

/* Counts the bits which are mismatched in 2 double values from double buffers
 *    - Load data into uint32_t to enable bitwise operations
 *    - Run XOR operation on the 2 uint32_t values for setting mismatch bits
 *    - Count the number of set bits in the XOR result to get the number of mismatch
 *  */
#define bitwise_compare_double(indexA, indexB)                     \
    {                                                              \
        uint32_t a[2], b[2];                                       \
        memcpy(&a, &((double *)A)[indexA], sizeof(double));        \
        memcpy(&b, &((double *)B)[indexB], sizeof(double));        \
        mismatchBits = mismatchBits + count_set_bits(a[0] ^ b[0]); \
        mismatchBits = mismatchBits + count_set_bits(a[1] ^ b[1]); \
    }

/* Counts the bits which are mismatched in 2 complex values from complex buffers
 *    - Load data into uint32_t to enable bitwise operations
 *    - Run XOR operation on the 2 uint32_t values for setting mismatch bits
 *    - Count the number of set bits in the XOR result to get the number of mismatch
 *  */
#define bitwise_compare_scomplex(indexA, indexB)                  \
    {                                                             \
        uint32_t a, b, c, d;                                      \
        memcpy(&a, &((scomplex *)A)[indexA].real, sizeof(float)); \
        memcpy(&b, &((scomplex *)B)[indexB].real, sizeof(float)); \
        memcpy(&c, &((scomplex *)A)[indexA].imag, sizeof(float)); \
        memcpy(&d, &((scomplex *)B)[indexB].imag, sizeof(float)); \
        mismatchBits = mismatchBits + count_set_bits(a ^ b);      \
        mismatchBits = mismatchBits + count_set_bits(c ^ d);      \
    }

/* Counts the bits which are mismatched in 2 double complex values from
 double complex buffers
 *    - Load data into uint32_t to enable bitwise operations
 *    - Run XOR operation on the 2 uint32_t values for setting mismatch bits
 *    - Count the number of set bits in the XOR result to get the number of mismatch
 *  */
#define bitwise_compare_dcomplex(indexA, indexB)                   \
    {                                                              \
        uint32_t a[2], b[2];                                       \
        memcpy(&a, &((dcomplex *)A)[indexA].real, sizeof(double)); \
        memcpy(&b, &((dcomplex *)B)[indexB].real, sizeof(double)); \
        mismatchBits = mismatchBits + count_set_bits(a[0] ^ b[0]); \
        mismatchBits = mismatchBits + count_set_bits(a[1] ^ b[1]); \
        memcpy(&a, &((dcomplex *)A)[indexA].imag, sizeof(double)); \
        memcpy(&b, &((dcomplex *)B)[indexB].imag, sizeof(double)); \
        mismatchBits = mismatchBits + count_set_bits(a[0] ^ b[0]); \
        mismatchBits = mismatchBits + count_set_bits(a[1] ^ b[1]); \
    }

/* Returns number of bitwise mismatch between 2 vectors */
int bitwise_compare_vector(integer datatype, integer m, void *A, void *B)
{
    int mismatchBits = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            /* Iterate through each element in the vector*/
            for(size_t i = 0; i < m; i++)
            {
                /* Compare the elements and add to mismatch bits if any */
                if(((float *)A)[i] != ((float *)B)[i])
                    bitwise_compare_float(i, i);
            }
            break;
        }
        case DOUBLE:
        {
            /* Iterate through each element in the vector */
            for(size_t i = 0; i < m; i++)
            {
                /* Compare the elements and add to mismatch bits if any */
                if(((double *)A)[i] != ((double *)B)[i])
                    bitwise_compare_double(i, i);
            }
            break;
        }
        case COMPLEX:
        {
            /* Iterate through each element in the vector */
            for(size_t i = 0; i < m; i++)
            {
                /* Compare the elements and add to mismatch bits if any */
                if(((scomplex *)A)[i].real != ((scomplex *)B)[i].real
                   || ((scomplex *)A)[i].imag != ((scomplex *)B)[i].imag)
                    bitwise_compare_scomplex(i, i);
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Iterate through each element in the vector */
            for(size_t i = 0; i < m; i++)
            {
                /* Compare the elements and add to mismatch bits if any */
                if(((dcomplex *)A)[i].real != ((dcomplex *)B)[i].real
                   || ((dcomplex *)A)[i].imag != ((dcomplex *)B)[i].imag)
                    bitwise_compare_dcomplex(i, i);
            }
            break;
        }
    }
    return mismatchBits;
}

/* Returns number of bitwise mismatch between 2 matrices */
int bitwise_compare_matrix(integer datatype, integer m, integer n, void *A, integer lda, void *B,
                           integer ldb)
{
    int mismatchBits = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            /* Iterate through each element in the matrix*/
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    /* Compare the elements and add to mismatch bits if any */
                    if(((float *)A)[i * lda + j] != ((float *)B)[i * ldb + j])
                        bitwise_compare_float(i * lda + j, i * ldb + j);
                }
            }
            break;
        }
        case DOUBLE:
        {
            /* Iterate through each element in the matrix */
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    /* Compare the elements and add to mismatch bits if any */
                    if(((double *)A)[i * lda + j] != ((double *)B)[i * ldb + j])
                        bitwise_compare_double(i * lda + j, i * ldb + j);
                }
            }
            break;
        }
        case COMPLEX:
        {
            /* Iterate through each element in the matrix */
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    /* Compare the elements and add to mismatch bits if any */
                    if(((scomplex *)A)[i * lda + j].real != ((scomplex *)B)[i * ldb + j].real
                       || ((scomplex *)A)[i * lda + j].imag != ((scomplex *)B)[i * ldb + j].imag)
                        bitwise_compare_scomplex(i * lda + j, i * ldb + j);
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Iterate through each element in the matrix */
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    /* Compare the elements and add to mismatch bits if any */
                    if(((dcomplex *)A)[i * lda + j].real != ((dcomplex *)B)[i * ldb + j].real
                       || ((dcomplex *)A)[i * lda + j].imag != ((dcomplex *)B)[i * ldb + j].imag)
                        bitwise_compare_dcomplex(i * lda + j, i * ldb + j);
                }
            }
            break;
        }
    }
    return mismatchBits;
}

/* Returns number of bitwise mismatch between 2 matrices ignoring the first nb diagonal elements*/
int bitwise_compare_matrix_no_nb_diag(integer datatype, integer m, integer n, integer nb, void *A,
                                      integer lda, void *B, integer ldb)
{
    int mismatchBits = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            /* Iterate through each element in the matrix */
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    /* Skip the first nb diagonal elements */
                    if(i == j && i < nb)
                    {
                        continue;
                    }
                    /* Compare the elements and add to mismatch bits if any */
                    if(((float *)A)[i * lda + j] != ((float *)B)[i * ldb + j])
                        bitwise_compare_float(i * lda + j, i * ldb + j);
                }
            }
            break;
        }
        case DOUBLE:
        {
            /* Iterate through each element in the matrix */
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    /* Skip the first nb diagonal elements */
                    if(i == j && i < nb)
                    {
                        continue;
                    }
                    /* Compare the elements and add to mismatch bits if any */
                    if(((double *)A)[i * lda + j] != ((double *)B)[i * ldb + j])
                        bitwise_compare_double(i * lda + j, i * ldb + j);
                }
            }
            break;
        }
        case COMPLEX:
        {
            /* Iterate through each element in the matrix */
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    /* Skip the first nb diagonal elements */
                    if(i == j && i < nb)
                    {
                        continue;
                    }
                    /* Compare the elements and add to mismatch bits if any */
                    if(((scomplex *)A)[i * lda + j].real != ((scomplex *)B)[i * ldb + j].real
                       || ((scomplex *)A)[i * lda + j].imag != ((scomplex *)B)[i * ldb + j].imag)
                        bitwise_compare_scomplex(i * lda + j, i * ldb + j);
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Iterate through each element in the matrix */
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    /* Skip the first nb diagonal elements */
                    if(i == j && i < nb)
                    {
                        continue;
                    }
                    /* Compare the elements and add to mismatch bits if any */
                    if(((dcomplex *)A)[i * lda + j].real != ((dcomplex *)B)[i * ldb + j].real
                       || ((dcomplex *)A)[i * lda + j].imag != ((dcomplex *)B)[i * ldb + j].imag)
                        bitwise_compare_dcomplex(i * lda + j, i * ldb + j);
                }
            }
            break;
        }
    }
    return mismatchBits;
}

/* Counting set bits using the canonical Kernighan’s bit counting or __builtin_popcount*/
int count_set_bits(uint32_t v)
{
#if defined(__GNUC__)
    return __builtin_popcount(v);
#else
    // Fallback to portable implementation
    int c = 0;
    while(v)
    {
        c += v & 1;
        v >>= 1;
    }
    return c;
#endif
}

/* Macro for binary division, this adds the MSB of a uint_32 to the right of the remainder and runs
 * XOR operation with the polynomial divisor followed by a left shift and does this for each bit in
 * the value to be divided*/
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
    switch(datatype)
    {
        case INTEGER:
        {
            /* Load data into variable of type uint32_t to enable bitwise operations */
            uint32_t temp;
            memcpy(&temp, &((integer *)A)[0], sizeof(integer));
            /* The first remainder is generated by using bitwise XOR operation with the polynomial
             * divisor*/
            if((temp & 0x80000000) == 0x80000000)
                remainder = temp ^ polynomial;
            else
                remainder = temp;
            remainder = remainder << 1;
            /* For each subsequent bit,
             * - Left shift the remainder to bring space for the new bit
             * - Add the new bit from next float to the remainder
             * - Left shift the next float and iteratively add the last bit for binary division*/
            for(size_t i = 1; i < m; i++)
            {
                // Pack bits from collected datatype to uint32_t for enabling binary addition
                uint32_t a;
                memcpy(&a, &((integer *)A)[i], sizeof(integer));

                BINARY_DIVISION_UINT32(a);
            }
            break;
        }
        case FLOAT:
        {
            /* Load data into variable of type uint32_t to enable bitwise operations */
            uint32_t temp;
            memcpy(&temp, &((float *)A)[0], sizeof(float));
            /* The first remainder is generated by using bitwise XOR operation with the polynomial
             * divisor*/
            if((temp & 0x80000000) == 0x80000000)
                remainder = temp ^ polynomial;
            else
                remainder = temp;
            remainder = remainder << 1;
            /* For each subsequent bit,
             * - Left shift the remainder to bring space for the new bit
             * - Add the new bit from next float to the remainder
             * - Left shift the next float and iteratively add the last bit for binary division*/
            for(size_t i = 1; i < m; i++)
            {
                // Pack bits from collected datatype to uint32_t for enabling binary addition
                uint32_t a;
                memcpy(&a, &((float *)A)[i], sizeof(float));

                BINARY_DIVISION_UINT32(a);
            }
            break;
        }
        case DOUBLE:
        {
            /* Load data into variable of type uint32_t to enable bitwise operations */
            uint32_t temp[2];
            memcpy(&temp, &((double *)A)[0], sizeof(double));
            /* The first remainder is generated by using bitwise XOR operation with the polynomial
             * divisor*/
            if((temp[0] & 0x80000000) == 0x80000000)
                remainder = temp[0] ^ polynomial;
            else
                remainder = temp[0];
            remainder = remainder << 1;
            /* For each subsequent bit,
             * - Left shift the remainder to bring space for the new bit
             * - Add the new bit from next float to the remainder
             * - Left shift the next float and iteratively add the last bit for binary division*/
            BINARY_DIVISION_UINT32(temp[1]);

            for(size_t i = 1; i < m; i++)
            {
                // Pack bits from collected datatype to uint32_t for enabling binary addition
                uint32_t a[2];
                memcpy(&a, &((double *)A)[i], sizeof(double));

                BINARY_DIVISION_UINT32(a[0]);
                BINARY_DIVISION_UINT32(a[1]);
            }
            break;
        }
        case COMPLEX:
        {
            /* Load data into variable of type uint32_t to enable bitwise operations */
            uint32_t tempa, tempb;
            memcpy(&tempa, &((scomplex *)A)[0].real, sizeof(float));
            memcpy(&tempb, &((scomplex *)A)[0].imag, sizeof(float));
            /* The first remainder is generated by using bitwise XOR operation with the polynomial
             * divisor*/
            if((tempa & 0x80000000) == 0x80000000)
                remainder = tempa ^ polynomial;
            else
                remainder = tempa;
            remainder = remainder << 1;
            /* For each subsequent bit,
             * - Left shift the remainder to bring space for the new bit
             * - Add the new bit from next float to the remainder
             * - Left shift the next float and iteratively add the last bit for binary division*/
            BINARY_DIVISION_UINT32(tempb);

            for(size_t i = 1; i < m; i++)
            {
                // Pack bits from collected datatype to uint32_t for enabling binary addition
                uint32_t a, b;
                memcpy(&a, &((scomplex *)A)[i].real, sizeof(float));
                memcpy(&b, &((scomplex *)A)[i].imag, sizeof(float));

                BINARY_DIVISION_UINT32(a);
                BINARY_DIVISION_UINT32(b);
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Load data into variable of type uint32_t to enable bitwise operations */
            uint32_t tempa[2], tempb[2];
            memcpy(&tempa, &((dcomplex *)A)[0].real, sizeof(double));
            memcpy(&tempb, &((dcomplex *)A)[0].imag, sizeof(double));
            /* The first remainder is generated by using bitwise XOR operation with the polynomial
             * divisor*/
            if((tempa[0] & 0x80000000) == 0x80000000)
                remainder = tempa[0] ^ polynomial;
            else
                remainder = tempa[0];
            remainder = remainder << 1;
            /* For each subsequent bit,
             * - Left shift the remainder to bring space for the new bit
             * - Add the new bit from next float to the remainder
             * - Left shift the next float and iteratively add the last bit for binary division*/
            BINARY_DIVISION_UINT32(tempa[1]);
            BINARY_DIVISION_UINT32(tempb[0]);
            BINARY_DIVISION_UINT32(tempb[1]);
            for(size_t i = 1; i < m; i++)
            {
                // Pack bits from collected datatype to uint32_t for enabling bitwise operations
                uint32_t a[2], b[2];
                memcpy(&a, &((dcomplex *)A)[i].real, sizeof(double));
                memcpy(&b, &((dcomplex *)A)[i].imag, sizeof(double));

                BINARY_DIVISION_UINT32(a[0]);
                BINARY_DIVISION_UINT32(a[1]);
                BINARY_DIVISION_UINT32(b[0]);
                BINARY_DIVISION_UINT32(b[1]);
            }
            break;
        }
    }
    return remainder;
}

/* The Cyclic redundancy check, checksum generation follows a bit by bit modulo-2 binary division
 * algorithm. This function returns the checksum generated for a matrix */
uint32_t generate_crc_matrix(integer datatype, integer m, integer n, void *A, integer lda)
{
    uint32_t remainder = 0;
    uint32_t polynomial = CRC32_POLYNOMIAL;
    switch(datatype)
    {
        case FLOAT:
        {
            /* Load data into variable of type uint32_t to enable bitwise operations */
            uint32_t temp;
            memcpy(&temp, &((float *)A)[0], sizeof(float));
            /* The first remainder is generated by using bitwise XOR operation with the polynomial
             * divisor*/
            if((temp & 0x80000000) == 0x80000000)
                remainder = temp ^ polynomial;
            else
                remainder = temp;
            remainder = remainder << 1;
            /* For each subsequent bit,
             * - Left shift the remainder to bring space for the new bit
             * - Add the new bit from next float to the remainder
             * - Left shift the next float and iteratively add the last bit for binary division*/
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    if((i * lda + j) != 0)
                    {
                        // Pack bits from collected datatype to uint32_t for enabling bitwise
                        // operations
                        uint32_t a;
                        memcpy(&a, &((float *)A)[i * lda + j], sizeof(float));
                        BINARY_DIVISION_UINT32(a);
                    }
                }
            }
            break;
        }
        case DOUBLE:
        {
            /* Load data into variable of type uint32_t to enable bitwise operations */
            uint32_t temp[2];
            memcpy(&temp, &((double *)A)[0], sizeof(double));
            /* The first remainder is generated by using bitwise XOR operation with the polynomial
             * divisor*/
            if((temp[0] & 0x80000000) == 0x80000000)
                remainder = temp[0] ^ polynomial;
            else
                remainder = temp[0];
            remainder = remainder << 1;
            BINARY_DIVISION_UINT32(temp[1]);
            /* For each subsequent bit,
             * - Left shift the remainder to bring space for the new bit
             * - Add the new bit from next float to the remainder
             * - Left shift the next float and iteratively add the last bit for binary division*/

            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    if((i * lda + j) != 0)
                    {
                        // Pack bits from collected datatype to uint32_t for enabling bitwise
                        // operations
                        uint32_t a[2];
                        memcpy(&a, &((double *)A)[i * lda + j], sizeof(double));

                        BINARY_DIVISION_UINT32(a[0]);
                        BINARY_DIVISION_UINT32(a[1]);
                    }
                }
            }
            break;
        }
        case COMPLEX:
        {
            /* Load data into variable of type uint32_t to enable bitwise operations */
            uint32_t tempa, tempb;
            memcpy(&tempa, &((scomplex *)A)[0].real, sizeof(float));
            memcpy(&tempb, &((scomplex *)A)[0].imag, sizeof(float));
            /* The first remainder is generated by using bitwise XOR operation with the polynomial
             * divisor*/
            if((tempa & 0x80000000) == 0x80000000)
                remainder = tempa ^ polynomial;
            else
                remainder = tempa;
            remainder = remainder << 1;
            BINARY_DIVISION_UINT32(tempb);
            /* For each subsequent bit,
             * - Left shift the remainder to bring space for the new bit
             * - Add the new bit from next float to the remainder
             * - Left shift the next float and iteratively add the last bit for binary division*/
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    if((i * lda + j) != 0)
                    {
                        uint32_t a, b;
                        memcpy(&a, &((scomplex *)A)[i * lda + j].real, sizeof(float));
                        memcpy(&b, &((scomplex *)A)[i * lda + j].imag, sizeof(float));

                        BINARY_DIVISION_UINT32(a);
                        BINARY_DIVISION_UINT32(b);
                    }
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Load data into variable of type uint32_t to enable bitwise operations */
            uint32_t tempa[2], tempb[2];
            memcpy(&tempa, &((dcomplex *)A)[0].real, sizeof(double));
            memcpy(&tempb, &((dcomplex *)A)[0].imag, sizeof(double));
            /* The first remainder is generated by using bitwise XOR operation with the polynomial
             * divisor*/
            if((tempa[0] & 0x80000000) == 0x80000000)
                remainder = tempa[0] ^ polynomial;
            else
                remainder = tempa[0];
            remainder = remainder << 1;
            BINARY_DIVISION_UINT32(tempa[1]);
            BINARY_DIVISION_UINT32(tempb[0]);
            BINARY_DIVISION_UINT32(tempb[1]);

            /* For each subsequent bit,
             * - Left shift the remainder to bring space for the new bit
             * - Add the new bit from next float to the remainder
             * - Left shift the next float and iteratively add the last bit for binary division*/
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    if((i * lda + j) != 0)
                    {
                        // Pack bits from collected datatype to uint32_t for enabling bitwise
                        // operations
                        uint32_t a[2], b[2];
                        memcpy(&b, &((dcomplex *)A)[i * lda + j].imag, sizeof(double));
                        memcpy(&a, &((dcomplex *)A)[i * lda + j].real, sizeof(double));

                        BINARY_DIVISION_UINT32(a[0]);
                        BINARY_DIVISION_UINT32(a[1]);
                        BINARY_DIVISION_UINT32(b[0]);
                        BINARY_DIVISION_UINT32(b[1]);
                    }
                }
            }
            break;
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
    switch(datatype)
    {
        case FLOAT:
        {
            /* Load data into variable of type uint32_t to enable bitwise operations */
            uint32_t temp;
            memcpy(&temp, &((float *)A)[1], sizeof(float));
            /* The first remainder is generated by using bitwise XOR operation with the polynomial
             * divisor*/
            if((temp & 0x80000000) == 0x80000000)
                remainder = temp ^ polynomial;
            else
                remainder = temp;
            remainder = remainder << 1;
            /* For each subsequent bit,
             * - Left shift the remainder to bring space for the new bit
             * - Add the new bit from next float to the remainder
             * - Left shift the next float and iteratively add the last bit for binary division*/
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    /* The first nb diagonals are not taken for consideration */
                    if(i == j && i < nb)
                    {
                        continue;
                    }
                    if((i * lda + j) != 1)
                    {
                        // Pack bits from collected datatype to uint32_t for enabling bitwise
                        // operations
                        uint32_t a;
                        memcpy(&a, &((float *)A)[i * lda + j], sizeof(float));

                        BINARY_DIVISION_UINT32(a);
                    }
                }
            }
            break;
        }
        case DOUBLE:
        {
            /* Load data into variable of type uint32_t to enable bitwise operations */
            uint32_t temp[2];
            memcpy(&temp, &((double *)A)[1], sizeof(double));
            /* The first remainder is generated by using bitwise XOR operation with the polynomial
             * divisor*/
            if((temp[0] & 0x80000000) == 0x80000000)
                remainder = temp[0] ^ polynomial;
            else
                remainder = temp[0];
            remainder = remainder << 1;
            BINARY_DIVISION_UINT32(temp[1]);

            /* For each subsequent bit,
             * - Left shift the remainder to bring space for the new bit
             * - Add the new bit from next float to the remainder
             * - Left shift the next float and iteratively add the last bit for binary division*/
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    /* The first nb diagonals are not taken for consideration */
                    if(i == j && i < nb)
                    {
                        continue;
                    }
                    if((i * lda + j) != 1)
                    {
                        // Pack bits from collected datatype to uint32_t for enabling bitwise
                        // operations
                        uint32_t a[2];
                        memcpy(&a, &((double *)A)[i * lda + j], sizeof(double));

                        BINARY_DIVISION_UINT32(a[0]);
                        BINARY_DIVISION_UINT32(a[1]);
                    }
                }
            }
            break;
        }
        case COMPLEX:
        {
            /* Load data into variable of type uint32_t to enable bitwise operations */
            uint32_t tempa, tempb;
            memcpy(&tempa, &((scomplex *)A)[1].real, sizeof(float));
            memcpy(&tempb, &((scomplex *)A)[1].imag, sizeof(float));
            /* The first remainder is generated by using bitwise XOR operation with the polynomial
             * divisor*/
            if((tempa & 0x80000000) == 0x80000000)
                remainder = tempa ^ polynomial;
            else
                remainder = tempa;
            remainder = remainder << 1;
            BINARY_DIVISION_UINT32(tempb);

            /* For each subsequent bit,
             * - Left shift the remainder to bring space for the new bit
             * - Add the new bit from next float to the remainder
             * - Left shift the next float and iteratively add the last bit for binary division*/
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    /* The first nb diagonals are not taken for consideration */
                    if(i == j && i < nb)
                    {
                        continue;
                    }
                    if((i * lda + j) != 1)
                    {
                        uint32_t a, b;
                        memcpy(&a, &((scomplex *)A)[i * lda + j].real, sizeof(float));
                        memcpy(&b, &((scomplex *)A)[i * lda + j].imag, sizeof(float));

                        BINARY_DIVISION_UINT32(a);
                        BINARY_DIVISION_UINT32(b);
                    }
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Load data into variable of type uint32_t to enable bitwise operations */
            uint32_t tempa[2], tempb[2];
            memcpy(&tempa, &((dcomplex *)A)[1].real, sizeof(double));
            memcpy(&tempb, &((dcomplex *)A)[1].imag, sizeof(double));
            /* The first remainder is generated by using bitwise XOR operation with the polynomial
             * divisor*/
            if((tempa[0] & 0x80000000) == 0x80000000)
                remainder = tempa[0] ^ polynomial;
            else
                remainder = tempa[0];
            remainder = remainder << 1;
            BINARY_DIVISION_UINT32(tempa[1]);
            BINARY_DIVISION_UINT32(tempb[0]);
            BINARY_DIVISION_UINT32(tempb[1]);

            /* For each subsequent bit,
             * - Left shift the remainder to bring space for the new bit
             * - Add the new bit from next float to the remainder
             * - Left shift the next float and iteratively add the last bit for binary division*/
            for(size_t i = 0; i < n; i++)
            {
                for(size_t j = 0; j < m; j++)
                {
                    /* The first nb diagonals are not taken for consideration */
                    if(i == j && i < nb)
                    {
                        continue;
                    }
                    if((i * lda + j) != 1)
                    {
                        // Pack bits from collected datatype to uint32_t for enabling bitwise
                        // operations
                        uint32_t a[2], b[2];
                        memcpy(&b, &((dcomplex *)A)[i * lda + j].imag, sizeof(double));
                        memcpy(&a, &((dcomplex *)A)[i * lda + j].real, sizeof(double));

                        BINARY_DIVISION_UINT32(a[0]);
                        BINARY_DIVISION_UINT32(a[1]);
                        BINARY_DIVISION_UINT32(b[0]);
                        BINARY_DIVISION_UINT32(b[1]);
                    }
                }
            }
            break;
        }
    }
    return remainder;
}