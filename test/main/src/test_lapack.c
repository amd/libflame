/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "ctype.h"
#include "test_lapack.h"
#include "test_routines.h"

// Global variables.
int n_threads = 1;

char *LINEAR_PARAMETERS_FILENAME = NULL;
char *SYM_EIG_PARAMETERS_FILENAME = NULL;
char *SVD_PARAMETERS_FILENAME = NULL;
char *NON_SYM_EIG_PARAMETERS_FILENAME = NULL;
char *AUX_PARAMETERS_FILENAME = NULL;

char fla_test_binary_name[MAX_BINARY_NAME_LENGTH + 1];
char fla_test_pass_string[MAX_PASS_STRING_LENGTH + 1];
char fla_test_warn_string[MAX_PASS_STRING_LENGTH + 1];
char fla_test_fail_string[MAX_PASS_STRING_LENGTH + 1];
char fla_test_invalid_string[MAX_PASS_STRING_LENGTH + 1];
char fla_test_storage_format_string[200];
char fla_test_stor_chars[NUM_STORAGE_CHARS + 1];
double ref_time_sec = 0.0;

integer g_total_tests;
integer g_total_failed_tests;
integer g_total_incomplete_tests;
integer g_tests_passed[4];
integer g_tests_failed[4];
integer g_tests_incomplete[4];

/* Flag to indicate lwork/liwork/lrwork availability status
 * <= 0 - To be calculated
 * > 0  - Use the value
 * */
integer g_lwork = -1;
integer g_liwork = -1;
integer g_lrwork = -1;
/* Variable to indicate the source of inputs
 * = 0 - Inputs are from command line
 * = 1 - Inputs are from config file
 * */
integer g_config_data = 0;
/* File pointer for external file which is used
 * to pass the input matrix values
 * */
FILE *g_ext_fptr = NULL;

/*
 * List of available statistics.
 * The order of the statistics in this list should match
 * the order of the fla_stat_type enum in test_lapack.h.
 */
fla_stat_info_t AVAILABLE_STATS[FLA_NUM_STATS]
    = {{FLA_MIN_STAT, "min", "MIN"},      {FLA_MAX_STAT, "max", "MAX"},
       {FLA_AVG_STAT, "avg", "AVG"},      {FLA_PERCENTILE_STAT, "p", "P"},
       {FLA_VARIANCE_STAT, "var", "VAR"}, {FLA_STDDEV_STAT, "stddev", "STD"}};

#define SKIP_EXTRA_LINE_READ() \
    eol = fgetc(fp);           \
    if(eol != '\n')            \
    fscanf(fp, "%*[^\n]\n")

#define CHECK_LINE_SKIP()                                           \
    eol = fgetc(fp);                                                \
    if((eol == '\r') || (eol == '\n'))                              \
    {                                                               \
        num_ranges = ((i + 1) < num_ranges) ? (i + 1) : num_ranges; \
        break;                                                      \
    }

int fla_check_cmd_config_dir(int argc, char **argv);

#if AOCL_FLA_SET_PROGRESS_ENABLE == 1
int aocl_fla_progress(const char *const api, const integer lenapi, const integer *const progress,
                      const integer *const current_thread, const integer *const total_threads)
{
    printf("In AOCL FLA  Progress thread  %" FT_IS ", at API  %s, progress  %" FT_IS
           " total threads= %" FT_IS "\n",
           *current_thread, api, *progress, *total_threads);
    return 0;
}
#endif

#if AOCL_FLA_SET_PROGRESS_ENABLE == 2
int test_progress(const char *const api, const integer lenapi, const integer *const progress,
                  const integer *const current_thread, const integer *const total_threads)
{
    printf("In AOCL Progress thread  %" FT_IS ", at API  %s, progress %" FT_IS
           " total threads= %" FT_IS " \n",
           *current_thread, api, *progress, *total_threads);
    return 0;
}
#endif

int main(int argc, char **argv)
{
    test_params_t params;
    integer arg_count = argc;
    bool status;

    /* Initialize some strings. */
    fla_test_init_strings();

    /* Check if multithread variable is set or not*/
    char *str = getenv("FLA_TEST_NUM_THREADS");

    if(str != NULL)
    {

        n_threads = atoi(str);

        if(n_threads < 1)
        {
            fla_test_output_error("Number of threads should be greater than or equal to 1 \n");
        }
    }

    /* Set default params */
    params.imatrix_char = '\0'; // Initialize imatrix_char to NULL
    params.interfacetype = LAPACK_TEST;
    params.matrix_major = LAPACK_COL_MAJOR;
    params.cli_print_header = 0;
    params.warmup_repeats = 0.;
    params.time_unit = FLA_TIME_UNIT_AUTO;
    params.outlier_multiplier = 0.0;
    params.dump_runtimes_file_name = NULL;

    status = fla_parse_cmdline_args(&arg_count, argv, &params);

    if(status == FALSE)
    {
        /* If parse error or --help then return */
        return -1;
    }

    /* Checking for the cmd option or config file option */
    int cmd_option = fla_check_cmd_config_dir(arg_count, argv);

    fla_test_runtime_ctx_init(&params);

    /* Check for Command line requests */
    if(cmd_option == 1)
    {
        g_config_data = 0;
        fla_test_execute_cli_api(arg_count, argv, &params);
    }
    else if(cmd_option == 0)
    {
#if ENABLE_AOCL_EXTENSION_APIS == 1
        printf(" AOCL-LAPACK version: %s\n", FLA_Get_AOCL_Version());
#endif
        g_config_data = 1;
        /* Copy the binary name to a global string so we can use it later. */
        strncpy(fla_test_binary_name, argv[0], MAX_BINARY_NAME_LENGTH);

        /* Read Linear API parameters from config file */
        fla_test_read_linear_param(LINEAR_PARAMETERS_FILENAME, &params);

        /* Read eigen parameters from config file */
        fla_test_read_sym_eig_params(SYM_EIG_PARAMETERS_FILENAME, &params);
        fla_test_read_non_sym_eig_params(NON_SYM_EIG_PARAMETERS_FILENAME, &params);

        /* Read SVD parameters from config file */
        fla_test_read_svd_params(SVD_PARAMETERS_FILENAME, &params);

        /* Read AUX parameters from config file */
        fla_test_read_aux_params(AUX_PARAMETERS_FILENAME, &params);

#if AOCL_FLA_SET_PROGRESS_ENABLE == 2
        aocl_fla_set_progress(test_progress);
#endif
        if((params.interfacetype == LAPACKE_ROW_TEST)
           || (params.interfacetype == LAPACKE_COLUMN_TEST))
        {
            fla_test_lapack_suite(LAPACKE_OPERATIONS_FILENAME, &params);
        }
        else
        {
            /* Test the LAPACK-level operations. */
            fla_test_lapack_suite(LAPACK_OPERATIONS_FILENAME, &params);
        }

        if(LINEAR_PARAMETERS_FILENAME)
            free(LINEAR_PARAMETERS_FILENAME);
        if(SYM_EIG_PARAMETERS_FILENAME)
            free(SYM_EIG_PARAMETERS_FILENAME);
        if(SVD_PARAMETERS_FILENAME)
            free(SVD_PARAMETERS_FILENAME);
        if(NON_SYM_EIG_PARAMETERS_FILENAME)
            free(NON_SYM_EIG_PARAMETERS_FILENAME);
        if(AUX_PARAMETERS_FILENAME)
            free(AUX_PARAMETERS_FILENAME);
    }
    else
    {
        return 0;
    }

    fla_test_runtime_ctx_free(&params);

    return 0;
}

/* Function to configure appropriate interface to test
   Returns true if interface is valid, returns false otherwise */
integer fla_check_interface(integer arg_count, char **argv, test_params_t *params)
{
    char *interface_test = "--interface=";
    char *row_major = "lapacke_row";
    char *column_major = "lapacke_column";
    char *cpp_test = "cpp";
    char *lapack_test = "lapack";
    char *interface_buff = NULL;
    int lapacke_major = LAPACK_COL_MAJOR;
    integer interfacetype = LAPACK_TEST;
    integer len_interface_test = strlen(interface_test);
    integer len_row_major = strlen(row_major);
    integer len_column_major = strlen(column_major);
    integer len_cpp_test = strlen(cpp_test);
    integer len_lapack_test = strlen(lapack_test);
    integer index;
    integer parse_status = 0;

    /* check all the input args excluding first argument test_lapack.x
       for '--interface=' string */
    for(index = 1; index < arg_count; index++)
    {
        if(!(strncmp(argv[index], interface_test, len_interface_test)))
        {
            parse_status = 1; /* argument is found */
            interface_buff = argv[index] + len_interface_test;

            for(int i = 0; i < strlen(argv[index]); i++)
            {
                interface_buff[i] = tolower(interface_buff[i]);
            }

            /* Check for specific interface like lapacke, cpp, lapack.*/
            if(!(strncmp(interface_buff, row_major, len_row_major))
               && (len_row_major == strlen(interface_buff)))
            {
                lapacke_major = LAPACK_ROW_MAJOR;
                interfacetype = LAPACKE_ROW_TEST;
            }
            else if(!(strncmp(interface_buff, column_major, len_column_major))
                    && (len_column_major == strlen(interface_buff)))
            {
                lapacke_major = LAPACK_COL_MAJOR;
                interfacetype = LAPACKE_COLUMN_TEST;
            }
            else if(!(strncmp(interface_buff, cpp_test, len_cpp_test))
                    && (len_cpp_test == strlen(interface_buff)))
            {
                interfacetype = LAPACK_CPP_TEST;
#if(!ENABLE_CPP_TEST)
                {
                    printf("\nError: ENABLE_CPP_TEST flag is disabled to use CPP interface, please "
                           "enable and rebuild.\n");
                    return -1; /* Invalid interface */
                }
#endif
            }
            else if(!(strncmp(interface_buff, lapack_test, len_lapack_test))
                    && (len_lapack_test == strlen(interface_buff)))
            /* assign default interface as lapack */
            {
                interfacetype = LAPACK_TEST;
            }
            else
            {
                printf("\nError: Interface '%s' is invalid,", interface_buff);
                printf(" Please provide valid interface: \n lapack, lapacke_row, lapacke_column, "
                       "cpp\n");
                return -1; /* Invalid interface */
            }
            break;
        }
    }
    /* Set these values only if interface paramter is found */
    if(parse_status)
    {
        params->interfacetype = interfacetype;
        params->matrix_major = lapacke_major;
    }
    return parse_status;
}

/* Function for checking cmd option or config file directory */
int fla_check_cmd_config_dir(int arg_count, char **argv)
{
    integer len_lin_file, len_eig_file, len_svd_file, len_eig_nsy_file, len_aux_file;
    int cmd_test_option = 0;
    char *config_dir = NULL;
    char *lin_file;
    char *eig_file;
    char *svd_file;
    char *eig_nsy_file;
    char *aux_file;
    char *config_opt = "--config-dir=";
    integer len_config_opt = strlen(config_opt);

    struct stat info;
    bool dir = 0;

    /*for default config*/
    if(arg_count == 1)
    {
        lin_file = "config/short/LIN_SLVR.dat";
        eig_file = "config/short/EIG_PARAMS.dat";
        svd_file = "config/short/SVD.dat";
        eig_nsy_file = "config/short/EIG_NSYM_PARAMS.dat";
        aux_file = "config/short/AUX_PARAMS.dat";

        len_lin_file = strlen(lin_file);
        len_eig_file = strlen(eig_file);
        len_svd_file = strlen(svd_file);
        len_eig_nsy_file = strlen(eig_nsy_file);
        len_aux_file = strlen(aux_file);

        LINEAR_PARAMETERS_FILENAME = (char *)malloc(len_lin_file + 1);
        SYM_EIG_PARAMETERS_FILENAME = (char *)malloc(len_eig_file + 1);
        SVD_PARAMETERS_FILENAME = (char *)malloc(len_svd_file + 1);
        NON_SYM_EIG_PARAMETERS_FILENAME = (char *)malloc(len_eig_nsy_file + 1);
        AUX_PARAMETERS_FILENAME = (char *)malloc(len_aux_file + 1);

        memcpy(LINEAR_PARAMETERS_FILENAME, lin_file, len_lin_file + 1);

        memcpy(SYM_EIG_PARAMETERS_FILENAME, eig_file, len_eig_file + 1);

        memcpy(SVD_PARAMETERS_FILENAME, svd_file, len_svd_file + 1);

        memcpy(NON_SYM_EIG_PARAMETERS_FILENAME, eig_nsy_file, len_eig_nsy_file + 1);

        memcpy(AUX_PARAMETERS_FILENAME, aux_file, len_aux_file + 1);

        return 0;
    }
    else if(arg_count == 2 && strlen(argv[1]) > len_config_opt)
    {
        /*checking config dir option or cmd*/
        if(!(strncmp(argv[1], config_opt, len_config_opt)))
        {
            config_dir = (char *)malloc(strlen(argv[1] + len_config_opt) + 1);
            memcpy(config_dir, argv[1] + len_config_opt, strlen(argv[1] + len_config_opt) + 1);

            char *c_str = "config/";
            lin_file = "/LIN_SLVR.dat";
            svd_file = "/SVD.dat";
            eig_file = "/EIG_PARAMS.dat";
            eig_nsy_file = "/EIG_NSYM_PARAMS.dat";
            aux_file = "/AUX_PARAMS.dat";

            integer len_cstr = strlen(c_str);
            integer len_dir = strlen(config_dir);
            len_lin_file = strlen(lin_file);
            len_eig_file = strlen(eig_file);
            len_svd_file = strlen(svd_file);
            len_eig_nsy_file = strlen(eig_nsy_file);
            len_aux_file = strlen(aux_file);

            char *dir_path = (char *)malloc(len_cstr + len_dir + 1);

            /*checking given Directory exist or not*/
            memcpy(dir_path, c_str, len_cstr);
            memcpy(dir_path + len_cstr, config_dir, len_dir + 1);

            if(stat(dir_path, &info) != 0)
            {
                printf("Error: '%s' directory not found under 'config' directory.  Exiting... \n",
                       config_dir);
                cmd_test_option = -1;
            }
            else if(info.st_mode & S_IFDIR) // S_ISDIR() doesn't exist on my windows
                dir = 1;
            else
            {
                printf("Error: '%s' directory not found under 'config' directory.  Exiting... \n",
                       config_dir);
                cmd_test_option = -1;
            }

            /*Reading the config directory*/
            if(dir)
            {
                LINEAR_PARAMETERS_FILENAME = (char *)malloc(len_cstr + len_dir + len_lin_file + 1);
                SYM_EIG_PARAMETERS_FILENAME = (char *)malloc(len_cstr + len_dir + len_eig_file + 1);
                SVD_PARAMETERS_FILENAME = (char *)malloc(len_cstr + len_dir + len_svd_file + 1);
                NON_SYM_EIG_PARAMETERS_FILENAME
                    = (char *)malloc(len_cstr + len_dir + len_eig_nsy_file + 1);
                AUX_PARAMETERS_FILENAME = (char *)malloc(len_cstr + len_dir + len_aux_file + 1);

                memcpy(LINEAR_PARAMETERS_FILENAME, c_str, len_cstr);
                memcpy(LINEAR_PARAMETERS_FILENAME + len_cstr, config_dir, len_dir);
                memcpy(LINEAR_PARAMETERS_FILENAME + len_cstr + len_dir, lin_file, len_lin_file + 1);

                memcpy(SYM_EIG_PARAMETERS_FILENAME, c_str, len_cstr);
                memcpy(SYM_EIG_PARAMETERS_FILENAME + len_cstr, config_dir, len_dir);
                memcpy(SYM_EIG_PARAMETERS_FILENAME + len_cstr + len_dir, eig_file,
                       len_eig_file + 1);

                memcpy(SVD_PARAMETERS_FILENAME, c_str, len_cstr);
                memcpy(SVD_PARAMETERS_FILENAME + len_cstr, config_dir, len_dir);
                memcpy(SVD_PARAMETERS_FILENAME + len_cstr + len_dir, svd_file, len_svd_file + 1);

                memcpy(NON_SYM_EIG_PARAMETERS_FILENAME, c_str, len_cstr);
                memcpy(NON_SYM_EIG_PARAMETERS_FILENAME + len_cstr, config_dir, len_dir);
                memcpy(NON_SYM_EIG_PARAMETERS_FILENAME + len_cstr + len_dir, eig_nsy_file,
                       len_eig_nsy_file + 1);

                memcpy(AUX_PARAMETERS_FILENAME, c_str, len_cstr);
                memcpy(AUX_PARAMETERS_FILENAME + len_cstr, config_dir, len_dir);
                memcpy(AUX_PARAMETERS_FILENAME + len_cstr + len_dir, aux_file, len_aux_file + 1);

                cmd_test_option = 0;
            }
            if(dir_path)
                free(dir_path);
        }
        else
        {
            cmd_test_option = 1;
        }
    }
    else
    {
        /*cmd option*/
        cmd_test_option = 1;
    }
    if(config_dir)
        free(config_dir);

    return cmd_test_option;
}

/* This function reads the operation file to execute selected LAPACK APIs*/
void fla_test_lapack_suite(char *input_filename, test_params_t *params)
{
    char buffer[INPUT_BUFFER_SIZE];
    integer check_flag;
    integer i, j, op, test_api_count, g_id, test_group_count;
    FILE *input_stream;

    g_total_tests = 0;
    g_total_failed_tests = 0;
    g_total_incomplete_tests = 0;
    g_tests_passed[0] = g_tests_passed[1] = g_tests_passed[2] = g_tests_passed[3] = 0;
    g_tests_failed[0] = g_tests_failed[1] = g_tests_failed[2] = g_tests_failed[3] = 0;
    g_tests_incomplete[0] = g_tests_incomplete[1] = g_tests_incomplete[2] = g_tests_incomplete[3]
        = 0;

    test_api_count = sizeof(API_test_functions) / sizeof(API_test_functions[0]);
    test_group_count = sizeof(API_test_group) / sizeof(API_test_group[0]);

    fla_test_output_info("\n");
    if((params->interfacetype == LAPACKE_ROW_TEST)
       || (params->interfacetype == LAPACKE_COLUMN_TEST))
    {
        fla_test_output_info("--- LAPACKE-level operation tests ---------------------\n");
        if(params->interfacetype == LAPACKE_ROW_TEST)
            fla_test_output_info("\nMatrix layout: Row-major\n");
        else
            fla_test_output_info("\nMatrix layout: Column-major\n");
    }
    else
    {
        fla_test_output_info("--- LAPACK-level operation tests ---------------------\n");
        if(params->interfacetype == LAPACK_CPP_TEST)
            fla_test_output_info("\nInterface: CPP\n");
    }
    fla_test_output_info("\n");

    // Attempt to open input file corresponding to input_filename as
    // read-only/binary.
    input_stream = fopen(input_filename, "rb");

    // Check for success.
    if(input_stream == NULL)
    {
        fla_test_output_error("Failed to open input file %s. Check existence and permissions.\n",
                              input_filename);
    }

    // Check for '2' option in input config and 3 for group testing
    check_flag = fla_test_check_run_only(input_stream, &op, buffer);

    if(check_flag == 3)
    {
        for(i = 0; i < test_group_count; i++)
        {
            fla_test_read_tests_for_op(input_stream, &op, buffer);
            if(op == 1)
            {
                g_id = fla_test_get_group_id(buffer);
                if(g_id >= 0)
                {
                    for(j = 0; j < test_api_count; j++)
                    {
                        if(API_test_functions[j].type == g_id)
                        {
                            API_test_functions[j].fp(1, NULL, params);
                        }
                    }
                }
            }
        }
    }
    else
    {
        while(fla_test_read_tests_for_op(input_stream, &op, buffer))
        {
            if(op == check_flag)
            {
                // Check if the specified API is supported in test suite
                for(i = 0; i < test_api_count; i++)
                {
                    if(!strcmp(API_test_functions[i].ops, buffer))
                    {
                        API_test_functions[i].fp(1, NULL, params);
                    }
                }

                op = 0;
            }
        }
    }
    fclose(input_stream);
    fla_test_print_summary();
}

void fla_test_output_op_struct(char *op_str, integer op)
{
    fla_test_output_info("%s LAPACK  %" FT_IS "\n", op_str, op);
}

/* THis function checks if operations file specified any run only option */
integer fla_test_check_run_only(FILE *input_stream, integer *op, char *buffer)
{
    integer run_only_flag = 1, i;
    integer test_group_count = sizeof(API_test_group) / sizeof(API_test_group[0]);

    /*checks if operations file specified any sub-group of API is enable*/
    for(i = 0; i < test_group_count; i++)
    {
        fla_test_read_tests_for_op(input_stream, op, buffer);
        if((*op) == 1)
        {
            run_only_flag = 3;
            fseek(input_stream, 0, SEEK_SET);
            return run_only_flag;
        }
    }

    while(fla_test_read_tests_for_op(input_stream, op, buffer))
    {
        run_only_flag = fla_max(run_only_flag, *op);
    }

    fseek(input_stream, 0, SEEK_SET);

    return run_only_flag;
}

/* This function reads group_id for given group name*/
integer fla_test_get_group_id(char *buffer)
{
    integer i;
    integer test_group_count = sizeof(API_test_group) / sizeof(API_test_group[0]);

    for(i = 0; i < test_group_count; i++)
    {
        if(!strcmp(buffer, API_test_group[i]))
            return i;
    }

    return -1;
}

/* This functiom extract enable option and API name from operation file */
integer fla_test_read_tests_for_op(FILE *input_stream, integer *op, char *buffer)
{
    char temp[INPUT_BUFFER_SIZE];

    // We want to read at least one line, so we use a do-while loop.
    do
    {
        // Read the next line into a temporary buffer and check success.
        if(fgets(temp, INPUT_BUFFER_SIZE - 1, input_stream) == NULL)
            return 0;
    }
    // We continue to read lines into buffer until the line is neither
    // commented nor blank.
    while(temp[0] == COMMENT_CHAR || temp[0] == '\n' || temp[0] == ' ' || temp[0] == '\t');

    // Save the string in temp, up to first white space character, into buffer.
    sscanf(temp, "%" FT_IS " %s", op, buffer);

    return 1;
}

/* This function reads operation file line by line*/
void fla_test_read_next_line(char *buffer, FILE *input_stream)
{
    char temp[INPUT_BUFFER_SIZE];

    // We want to read at least one line, so we use a do-while loop.
    do
    {
        // Read the next line into a temporary buffer and check success.
        if(fgets(temp, INPUT_BUFFER_SIZE - 1, input_stream) == NULL)
        {
            if(feof(input_stream))
                fla_test_output_error("Error reading input file: encountered unexpected EOF.");
            else
                fla_test_output_error("Error (non-EOF) reading input file.");
        }
    }
    // We continue to read lines into buffer until the line is neither
    // commented nor blank.
    while(temp[0] == COMMENT_CHAR || temp[0] == '\n' || temp[0] == ' ' || temp[0] == '\t');

    // Save the string in temp, up to first white space character, into buffer.
    sscanf(temp, "%s ", buffer);
}

#define READ_CONFIG_PARAM_INT(x)        \
    fscanf(fp, "%s", &line[0]);         \
    for(i = 0; i < NUM_SUB_TESTS; i++)  \
    {                                   \
        fscanf(fp, "%" FT_IS "", &(x)); \
        CHECK_LINE_SKIP();              \
    }

#define READ_CONFIG_PARAM_FLT(x)       \
    fscanf(fp, "%s", &line[0]);        \
    for(i = 0; i < NUM_SUB_TESTS; i++) \
    {                                  \
        fscanf(fp, "%f", &(x));        \
        CHECK_LINE_SKIP();             \
    }

#define READ_CONFIG_PARAM_DBL(x)       \
    fscanf(fp, "%s", &line[0]);        \
    for(i = 0; i < NUM_SUB_TESTS; i++) \
    {                                  \
        fscanf(fp, "%lf", &(x));       \
        CHECK_LINE_SKIP();             \
    }

#define READ_CONFIG_PARAM_STR(x)       \
    fscanf(fp, "%s", str);             \
    for(i = 0; i < NUM_SUB_TESTS; i++) \
    {                                  \
        fscanf(fp, "%s", str);         \
        (x) = *str;                    \
        CHECK_LINE_SKIP();             \
    }
#define PARSE_CONFIG_DATATYPES(x)                                          \
    fscanf(fp, "%s", &line[0]); /* num data types */                       \
    for(i = 0; i < NUM_SUB_TESTS; i++)                                     \
    {                                                                      \
        fscanf(fp, "%s", str); /* num data types */                        \
        for(j = 0; j < NUM_SUB_TESTS; j++)                                 \
        {                                                                  \
            x[j].data_types_char[i] = *str;                                \
            x[j].data_types[i] = get_datatype(*str);                       \
        }                                                                  \
        eol = fgetc(fp);                                                   \
        if((eol == '\r') || (eol == '\n'))                                 \
        {                                                                  \
            ndata_types = ((i + 1) < ndata_types) ? (i + 1) : ndata_types; \
            break;                                                         \
        }                                                                  \
    }                                                                      \
    for(i = 0; i < NUM_SUB_TESTS; i++)                                     \
    {                                                                      \
        x[i].num_data_types = ndata_types;                                 \
    }

/* This function reads parameters needed for Linear solver APIs
   from the config settings file 'LIN_SLVR.dat' and saves in the
   'lin_solver_paramslist' structure array   */
void fla_test_read_linear_param(const char *file_name, test_params_t *params)
{
    FILE *fp;
    integer i, j;
    char line[20];
    char *str;
    char eol;
    integer num_tests, ndata_types;
    integer num_ranges;

    str = &line[0];
    fp = fopen(file_name, "r");
    if(fp == NULL)
    {
        printf("Error: Lin solver config file missing. Exiting.. \n");
        exit(-1);
    }

    /* Read the number of Tests */
    fscanf(fp, "%s", &line[0]);
    fscanf(fp, "%" FT_IS "", &num_tests);

    num_ranges = num_tests;
    for(i = 0; i < NUM_SUB_TESTS; i++)
    {
        params->lin_solver_paramslist[i].num_tests = num_tests;
    }

    /* Range Start */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].m_range_start);
    /* Range End */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].m_range_end);
    /* Range_step_size */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].m_range_step_size);
    /* Range_start */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].n_range_start);
    /* Range_end */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].n_range_end);
    /* Range_step_size */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].n_range_step_size);
    /* leading dimension for A */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].lda);
    /* leading dimension for B */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].ldb);
    /* leading dimension for C */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].ldc);
    /* leading dimension for Q */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].ldq);
    /* leading dimension for Z */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].ldz);
    /* leading dimension LDAB */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].ldab);

    for(i = 0; i < NUM_SUB_TESTS; i++)
    {
        params->lin_solver_paramslist[i].num_ranges = num_ranges;
    }

    /* Number of repeats */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].num_repeats);

    /* datatypes specified */
    ndata_types = NUM_SUB_TESTS;
    PARSE_CONFIG_DATATYPES(params->lin_solver_paramslist);

    /* Matrix Layout (row or col major) */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].matrix_layout);
    /* trans */
    READ_CONFIG_PARAM_STR(params->lin_solver_paramslist[i].transr);
    /* uplo */
    READ_CONFIG_PARAM_STR(params->lin_solver_paramslist[i].Uplo);
    /* compq_gghrd */
    READ_CONFIG_PARAM_STR(params->lin_solver_paramslist[i].compq_gghrd);
    /* compz_gghrd */
    READ_CONFIG_PARAM_STR(params->lin_solver_paramslist[i].compz_gghrd);
    /* nrhs */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].nrhs);
    /* ncolm */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].ncolm);
    /* kl */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].kl);
    /* ku */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].ku);
    /* kd */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].kd);
    /* diag */
    READ_CONFIG_PARAM_STR(params->lin_solver_paramslist[i].diag);
    /* fact */
    READ_CONFIG_PARAM_STR(params->lin_solver_paramslist[i].fact);
    /* equed */
    READ_CONFIG_PARAM_STR(params->lin_solver_paramslist[i].equed);
    /* symm */
    READ_CONFIG_PARAM_STR(params->lin_solver_paramslist[i].symm);
    /* equed_porfsx */
    READ_CONFIG_PARAM_STR(params->lin_solver_paramslist[i].equed_porfsx);
    /* n_err_bnds_porfsx */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].n_err_bnds_porfsx);
    /* nparams_porfsx */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].nparams_porfsx);
    /* norm_gbcon */
    READ_CONFIG_PARAM_STR(params->lin_solver_paramslist[i].norm_gbcon);
    /* kl_gbcon */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].kl_gbcon);
    /* ku_gbcon */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].ku_gbcon);
    /* ldab_gbcon */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].ldab_gbcon);
    /* ilo */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].ilo);
    /* ihi */
    READ_CONFIG_PARAM_INT(params->lin_solver_paramslist[i].ihi);
    /* rcond */
    READ_CONFIG_PARAM_DBL(params->lin_solver_paramslist[i].rcond);
    /* side */
    READ_CONFIG_PARAM_STR(params->lin_solver_paramslist[i].side);
    /* solver_threshold */
    READ_CONFIG_PARAM_FLT(params->lin_solver_paramslist[i].solver_threshold);

    fclose(fp);
}

/* This function reads parameters needed for Eigen APIs
   from the config settings file 'EIG_PARAMS.dat' and saves in the
   'params->eig_sym_paramslist' structure array   */
void fla_test_read_sym_eig_params(const char *file_name, test_params_t *params)
{
    FILE *fp;
    integer i, j;
    char line[20], eol;
    char *str, c[20];
    integer num_tests;
    integer num_ranges;
    integer ndata_types = NUM_SUB_TESTS;

    str = &c[0];
    fp = fopen(file_name, "r");
    if(fp == NULL)
    {
        printf("Error: Symmetric EIG params config file missing. Exiting.. \n");
        exit(-1);
    }
    str = &line[0];

    /* Read the number of Ranges */
    fscanf(fp, "%s", &line[0]);
    fscanf(fp, "%" FT_IS "", &num_tests);
    for(i = 0; i < NUM_SUB_TESTS; i++)
    {
        params->eig_sym_paramslist[i].num_tests = num_tests;
    }

    num_ranges = num_tests;

    /* m_range_start */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].m_range_start);
    /* m_range_end */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].m_range_end);
    /* m_range_step_size */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].m_range_step_size);
    /* n_range_start */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].n_range_start);
    /* n_range_end */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].n_range_end);
    /* n_range_step_size */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].n_range_step_size);

    for(i = 0; i < NUM_SUB_TESTS; i++)
    {
        params->eig_sym_paramslist[i].num_ranges = num_ranges;
    }

    /* num_repeats */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].num_repeats);

    /* datatypes specified */
    ndata_types = NUM_SUB_TESTS;
    PARSE_CONFIG_DATATYPES(params->eig_sym_paramslist);

    /* matrix_layout */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].matrix_layout);
    /* trans */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].trans);
    /* uplo */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].uplo);
    /* job */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].job);
    /* jobz */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].jobz);
    /* job_seqr */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].job_seqr);
    /* vect */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].vect);
    /* nrhs */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].nrhs);
    /* lda */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].lda);
    /* ldb */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].ldb);
    /* ldz */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].ldz);
    /* ldq */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].ldq);
    /* nb */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].nb);
    /* ldt */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].ldt);
    /* k */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].k);
    /* isgn */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].isgn);
    /* compq_hgeqz */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].compq_hgeqz);
    /* compz_hgeqz  */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].compz_hgeqz);
    /* compz */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].compz);
    /* compz_hseqr */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].compz_hseqr);
    /* kb */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].kb);
    /* itype */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].itype);
    /* vect_rd */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].vect_rd);
    /* side */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].side);
    /* eigsrc */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].eigsrc);
    /* initv */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].initv);
    /* norm */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].norm);
    /* diag */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].diag);
    /* storev */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].storev);
    /* tsize */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].tsize);
    /* ilo */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].ilo);
    /* ihi */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].ihi);
    /* Range is used to select the range of eigen values to be generated */
    READ_CONFIG_PARAM_STR(params->eig_sym_paramslist[i].range_x);
    /* Index of the smallest eigen value to be returned */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].IL);
    /* Index of the largest eigen value to be returned */
    READ_CONFIG_PARAM_INT(params->eig_sym_paramslist[i].IU);
    /* Lower bound of the interval to be searched for eigen values */
    READ_CONFIG_PARAM_FLT(params->eig_sym_paramslist[i].VL);
    /* Upper bound of the interval to be searched for eigen values */
    READ_CONFIG_PARAM_FLT(params->eig_sym_paramslist[i].VU);
    /* The absolute error tolerance for the eigen values */
    READ_CONFIG_PARAM_FLT(params->eig_sym_paramslist[i].abstol);
    /* threshold_value */
    READ_CONFIG_PARAM_FLT(params->eig_sym_paramslist[i].threshold_value);

    fclose(fp);
}

/* This function reads parameters needed for Non symmetric Eigen APIs
   from the config settings file 'EIG_NSYM_PARAMS.dat' and saves in the
   'params->eig_non_sym_paramslist' structure array   */
/* This function reads parameters needed for Non symmetric Eigen APIs
   from the config settings file 'EIG_NSYM_PARAMS.dat' and saves in the
   'params->eig_non_sym_paramslist' structure array   */

void fla_test_read_non_sym_eig_params(const char *file_name, test_params_t *params)
{
    FILE *fp;
    integer i, j;
    char line[20], eol;
    char *str = &line[0];
    integer num_tests;
    integer num_ranges;
    integer ndata_types = NUM_SUB_TESTS;
    fp = fopen(file_name, "r");
    if(fp == NULL)
    {
        printf("Error: EIG non symmetric API params config file missing. Exiting.. \n");
        exit(-1);
    }

    /* Read the number of Tests */
    fscanf(fp, "%s", &line[0]);
    fscanf(fp, "%" FT_IS "", &num_tests);
    for(i = 0; i < NUM_SUB_TESTS; i++)
    {
        params->eig_non_sym_paramslist[i].num_tests = num_tests;
    }

    num_ranges = num_tests;

    /* m_range_start */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].m_range_start);
    /* m_range_end */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].m_range_end);
    /* m_range_step_size */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].m_range_step_size);
    /* n_range_start */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].n_range_start);
    /* n_range_end */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].n_range_end);
    /* n_range_step_size */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].n_range_step_size);

    for(i = 0; i < NUM_SUB_TESTS; i++)
    {
        params->eig_non_sym_paramslist[i].num_ranges = num_ranges;
    }

    /* lda */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].lda);
    /* ldb */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].ldb);
    /* ldvl */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].ldvl);
    /* ldvr */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].ldvr);
    /* num_repeats */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].num_repeats);

    /* datatypes specified */
    ndata_types = NUM_SUB_TESTS;
    PARSE_CONFIG_DATATYPES(params->eig_non_sym_paramslist);

    /* matrix_layout */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].matrix_layout);
    /* howmny */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].howmny);
    /* initv */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].initv);
    /* job_seqr */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].job_seqr);
    /* eigsrc */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].eigsrc);
    /* initv */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].initv);
    /* job */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].job);
    /* howmny_trsna */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].howmny_trsna);
    /* job_trsen */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].job_trsen);
    /* compq */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].compq);
    /* Reading config params for 'trsyl' API  */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].trana_real);
    /* trana_complex */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].trana_complex);
    /* tranb_real */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].tranb_real);
    /* tranb_complex */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].tranb_complex);
    /* isgn */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].isgn);
    /* gghrd_threshold */
    READ_CONFIG_PARAM_FLT(params->eig_non_sym_paramslist[i].gghrd_threshold);
    /* ggbal_threshold */
    READ_CONFIG_PARAM_FLT(params->eig_non_sym_paramslist[i].ggbal_threshold);
    /* GenNonSymEigProblem_threshold */
    READ_CONFIG_PARAM_FLT(params->eig_non_sym_paramslist[i].GenNonSymEigProblem_threshold);
    /* side_tgevc */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].side_tgevc);
    /* jobvsl */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].jobvsl);
    /* jobvsr */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].jobvsr);
    /* sort */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].sort);
    /* sense_ggesx */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].sense_ggesx);
    /* balance_ggevx */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].balance_ggevx);
    /* sense_ggevx */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].sense_ggevx);
    /* sort_gees */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].sort_gees);
    /* wantz */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].wantz);
    /* wantq */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].wantq);
    /* tgsen_ijob */
    READ_CONFIG_PARAM_INT(params->eig_non_sym_paramslist[i].tgsen_ijob);
    /* unmhr_trans */
    READ_CONFIG_PARAM_STR(params->eig_non_sym_paramslist[i].unmhr_trans);

    fclose(fp);
}

/* This function reads parameters needed for SVD APIs
   from the config settings file 'SVD.dat' and saves in the
   'params->svd_paramslist' structure array   */
void fla_test_read_svd_params(const char *file_name, test_params_t *params)
{
    FILE *fp;
    integer i, j;
    char line[25];
    char *str;
    char eol;
    integer num_tests;
    integer ndata_types;
    integer num_ranges;

    str = &line[0];
    fp = fopen(file_name, "r");
    if(fp == NULL)
    {
        printf("Error: SVD config file missing. Exiting.. \n");
        exit(-1);
    }

    fscanf(fp, "%s", &line[0]);
    fscanf(fp, "%" FT_IS "", &num_tests);
    for(i = 0; i < NUM_SUB_TESTS; i++)
    {
        params->svd_paramslist[i].num_tests = num_tests;
    }

    num_ranges = num_tests;
    /* m_range_start */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].m_range_start);
    /* m_range_end */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].m_range_end);
    /* m_range_step_size */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].m_range_step_size);
    /* n_range_start */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].n_range_start);
    /* n_range_end */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].n_range_end);
    /* n_range_step_size */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].n_range_step_size);

    for(i = 0; i < NUM_SUB_TESTS; i++)
    {
        params->svd_paramslist[i].num_ranges = num_ranges;
    }

    /* lda */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].lda);
    /* ldu */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].ldu);
    /* ldvt */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].ldvt);
    /* num_repeats */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].num_repeats);

    /* datatypes specified */
    ndata_types = NUM_SUB_TESTS;
    PARSE_CONFIG_DATATYPES(params->svd_paramslist);

    /* matrix_layout */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].matrix_layout);
    /* jobu */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobu);
    /* jobv */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobv);
    /* jobq */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobq);
    /* m */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].m);
    /* p */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].p);
    /* n */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].n);
    /* tola */
    READ_CONFIG_PARAM_FLT(params->svd_paramslist[i].tola);
    /* tolb */
    READ_CONFIG_PARAM_FLT(params->svd_paramslist[i].tolb);
    /* jobu_gesvd */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobu_gesvd);
    /* jobvt_gesvd */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobvt_gesvd);
    /* joba_gejsv */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].joba_gejsv);
    /* jobu_gejsv */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobu_gejsv);
    /* jobv_gejsv */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobv_gejsv);
    /* jobr_gejsv */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobr_gejsv);
    /* jobt_gejsv */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobt_gejsv);
    /* jobp_gejsv */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobp_gejsv);
    /* m_gejsv */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].m_gejsv);
    /* n_gejsv */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].n_gejsv);
    /* joba_gesvj */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].joba_gesvj);
    /* jobu_gesvj */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobu_gesvj);
    /* jobv_gesvj */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobv_gesvj);
    /* m_gesvj */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].m_gesvj);
    /* n_gesvj */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].n_gesvj);
    /* mv_gesvj */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].mv_gesvj);
    /* ctol_gesvj */
    READ_CONFIG_PARAM_FLT(params->svd_paramslist[i].ctol_gesvj);
    /* jobu_gesvdx */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobu_gesvdx);
    /* jobvt_gesvdx */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobvt_gesvdx);
    /* range_gesvdx */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].range_gesvdx);
    /* il */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].il);
    /* iu */
    READ_CONFIG_PARAM_INT(params->svd_paramslist[i].iu);
    /* vl */
    READ_CONFIG_PARAM_FLT(params->svd_paramslist[i].vl);
    /* vu */
    READ_CONFIG_PARAM_FLT(params->svd_paramslist[i].vu);
    /* joba_gesvdq */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].joba_gesvdq);
    /* jobu_gesvdq */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobu_gesvdq);
    /* jobv_gesvdq */
    READ_CONFIG_PARAM_STR(params->svd_paramslist[i].jobv_gesvdq);
    /* svd_threshold */
    READ_CONFIG_PARAM_FLT(params->svd_paramslist[i].svd_threshold);

    fclose(fp);
}

/* This function reads parameters needed for aux APIs
   from the config settings file 'AUX_PARAMS.dat' and saves in the
   'params->aux_paramslist' structure array   */
void fla_test_read_aux_params(const char *file_name, test_params_t *params)
{
    FILE *fp;
    integer i, j;
    char line[25];
    char *str;
    char eol;
    integer num_tests;
    integer ndata_types;
    integer num_ranges;
    integer len;

    str = &line[0];
    fp = fopen(file_name, "r");
    if(fp == NULL)
    {
        printf("Error: aux config file missing. Exiting.. \n");
        exit(-1);
    }

    fscanf(fp, "%s", &line[0]);
    fscanf(fp, "%" FT_IS "", &num_tests);
    for(i = 0; i < NUM_SUB_TESTS; i++)
    {
        params->aux_paramslist[i].num_tests = num_tests;
    }

    num_ranges = num_tests;
    /* m_range_start */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].m_range_start);
    /* m_range_end */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].m_range_end);
    /* m_range_step_size */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].m_range_step_size);
    /* n_range_start */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].n_range_start);
    /* n_range_end */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].n_range_end);
    /* n_range_step_size */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].n_range_step_size);

    for(i = 0; i < NUM_SUB_TESTS; i++)
    {
        params->aux_paramslist[i].num_ranges = num_ranges;
    }

    /* nb */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].nb);
    /* lda */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].lda);
    /* ldx */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].ldx);
    /* ldy */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].ldy);
    /* incx */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].incx);
    /* incy */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].incy);
    /* incx_larfg */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].incx_larfg);
    /* side */
    READ_CONFIG_PARAM_STR(params->aux_paramslist[i].side);
    /* incv */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].incv);
    /* ldc */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].ldc);
    /* num_repeats */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].num_repeats);

    /* datatypes specified */
    ndata_types = NUM_SUB_TESTS;
    PARSE_CONFIG_DATATYPES(params->aux_paramslist);

    /* matrix_layout */
    READ_CONFIG_PARAM_INT(params->aux_paramslist[i].matrix_layout);
    /* aux_threshold */
    READ_CONFIG_PARAM_FLT(params->aux_paramslist[i].aux_threshold);

    fscanf(fp, "%s", &line[0]); // Norm types
    for(i = 0; i < NUM_SUB_TESTS; i++)
    {
        fscanf(fp, "%s", str);
        len = strlen(str);
        len = fla_min(len, MAX_NUM_NORMTYPES);
        for(j = 0; j < len; j++)
        {
            params->aux_paramslist[i].norm_types_str[j] = str[j];
        }
        for(j = len; j < MAX_NUM_NORMTYPES; j++)
        {
            params->aux_paramslist[i].norm_types_str[j] = '\0';
        }

        CHECK_LINE_SKIP();
    }

    fclose(fp);
}

void fla_test_execute_cli_api(integer argc, char **argv, test_params_t *params)
{
    integer i, test_api_count;
    char s_name[MAX_FUNC_STRING_LENGTH];

    FLA_PRINT_HEADER_CLI;

    test_api_count = sizeof(API_test_functions) / sizeof(API_test_functions[0]);

    /* Check if the specified API is supported in test suite */
    if(strlen(argv[1]) < MAX_FUNC_STRING_LENGTH)
    {
        strcpy(s_name, argv[1]);
        for(i = 0; i < strlen(s_name); i++)
        {
            s_name[i] = tolower(s_name[i]);
        }
        for(i = 0; i < test_api_count; i++)
        {
            if(!strcmp(API_test_functions[i].ops, s_name))
            {
                API_test_functions[i].fp(argc, argv, params);
                break;
            }
        }
    }

    if(g_total_tests == 0)
        printf("\nNo test was run, give valid arguments\n");
}

char *fla_test_get_string_for_result(double residual, integer datatype, double thresh)
{
    char *r_val;

    if(datatype == FLOAT)
    {
        if(residual == DBL_MIN)
            r_val = fla_test_invalid_string;
        else if((residual > thresh) || (isnan(residual)))
            r_val = fla_test_fail_string;
        else
            r_val = fla_test_pass_string;
    }
    else if(datatype == DOUBLE)
    {
        if(residual == DBL_MIN)
            r_val = fla_test_invalid_string;
        else if((residual > thresh) || (isnan(residual)))
            r_val = fla_test_fail_string;
        else
            r_val = fla_test_pass_string;
    }
    else if(datatype == COMPLEX)
    {
        if(residual == DBL_MIN)
            r_val = fla_test_invalid_string;
        else if((residual > thresh) || (isnan(residual)))
            r_val = fla_test_fail_string;
        else
            r_val = fla_test_pass_string;
    }
    else
    {
        if(residual == DBL_MIN)
            r_val = fla_test_invalid_string;
        else if((residual > thresh) || (isnan(residual)))
            r_val = fla_test_fail_string;
        else
            r_val = fla_test_pass_string;
    }

    return r_val;
}

void fla_test_init_strings(void)
{
    sprintf(fla_test_pass_string, "PASS");
    sprintf(fla_test_warn_string, "MARGINAL");
    sprintf(fla_test_fail_string, "FAIL");
    sprintf(fla_test_invalid_string, "INVALID_PARAM");
    sprintf(fla_test_storage_format_string,
            "Row(r) and General(g) storage format is not supported by External LAPACK interface");
    sprintf(fla_test_stor_chars, STORAGE_SCHEME_CHARS);
}

void fla_test_op_driver(char *func_str, integer sqr_inp, test_params_t *params, integer api_type,
                        void (*f_exp)(char *, // API_test string
                                      test_params_t *, // params
                                      integer, // datatype
                                      integer, // p_cur
                                      integer, // q_cur
                                      integer, // pci (param combo counter)
                                      integer, // n_repeats
                                      integer))
{
    integer n_datatypes = params->n_datatypes;
    integer n_repeats, ith;
    integer num_ranges, range_loop_counter;
    integer p_first, p_max, p_inc;
    integer q_first, q_max, q_inc;
    integer dt, p_cur, q_cur, einfo = 0;
    char datatype_char;
    integer datatype;
    double *perf = (double *)malloc(n_threads * sizeof(double));
    double *time = (double *)malloc(n_threads * sizeof(double));

    fla_test_print_header(params);

    switch(api_type)
    {
        case LIN:
            num_ranges = params->lin_solver_paramslist[0].num_ranges;
            break;

        case EIG_SYM:
            num_ranges = params->eig_sym_paramslist[0].num_ranges;
            break;

        case EIG_NSYM:
            num_ranges = params->eig_non_sym_paramslist[0].num_ranges;
            break;

        case SVD:
            num_ranges = params->svd_paramslist[0].num_ranges;
            break;

        case AUX:
            num_ranges = params->aux_paramslist[0].num_ranges;
            break;

        default:
            fla_test_output_error("Invalid API type. Exiting...\n");
            return;
    }

    // Loop over the requested ranges.
    for(range_loop_counter = 0; range_loop_counter < num_ranges; ++range_loop_counter)
    {
        switch(api_type)
        {
            case LIN:
                p_first = params->lin_solver_paramslist[range_loop_counter].m_range_start;
                p_max = params->lin_solver_paramslist[range_loop_counter].m_range_end;
                p_inc = params->lin_solver_paramslist[range_loop_counter].m_range_step_size;
                q_first = params->lin_solver_paramslist[range_loop_counter].n_range_start;
                q_max = params->lin_solver_paramslist[range_loop_counter].n_range_end;
                q_inc = params->lin_solver_paramslist[range_loop_counter].n_range_step_size;
                params->datatype = params->lin_solver_paramslist[range_loop_counter].data_types;
                params->datatype_char
                    = params->lin_solver_paramslist[range_loop_counter].data_types_char;
                n_repeats = params->lin_solver_paramslist[range_loop_counter].num_repeats;
                n_datatypes = params->lin_solver_paramslist[range_loop_counter].num_data_types;
                break;

            case EIG_SYM:
                p_first = params->eig_sym_paramslist[range_loop_counter].m_range_start;
                p_max = params->eig_sym_paramslist[range_loop_counter].m_range_end;
                p_inc = params->eig_sym_paramslist[range_loop_counter].m_range_step_size;
                q_first = p_first;
                q_max = p_max;
                q_inc = p_inc;
                params->datatype = params->eig_sym_paramslist[range_loop_counter].data_types;
                params->datatype_char
                    = params->eig_sym_paramslist[range_loop_counter].data_types_char;
                n_repeats = params->eig_sym_paramslist[range_loop_counter].num_repeats;
                n_datatypes = params->eig_sym_paramslist[range_loop_counter].num_data_types;
                break;

            case EIG_NSYM:
                p_first = params->eig_non_sym_paramslist[range_loop_counter].m_range_start;
                p_max = params->eig_non_sym_paramslist[range_loop_counter].m_range_end;
                p_inc = params->eig_non_sym_paramslist[range_loop_counter].m_range_step_size;
                q_first = params->eig_non_sym_paramslist[range_loop_counter].n_range_start;
                q_max = params->eig_non_sym_paramslist[range_loop_counter].n_range_end;
                q_inc = params->eig_non_sym_paramslist[range_loop_counter].n_range_step_size;
                params->datatype = params->eig_non_sym_paramslist[range_loop_counter].data_types;
                params->datatype_char
                    = params->eig_non_sym_paramslist[range_loop_counter].data_types_char;
                n_repeats = params->eig_non_sym_paramslist[range_loop_counter].num_repeats;
                n_datatypes = params->eig_non_sym_paramslist[range_loop_counter].num_data_types;
                break;

            case SVD:
                p_first = params->svd_paramslist[range_loop_counter].m_range_start;
                p_max = params->svd_paramslist[range_loop_counter].m_range_end;
                p_inc = params->svd_paramslist[range_loop_counter].m_range_step_size;
                q_first = params->svd_paramslist[range_loop_counter].n_range_start;
                q_max = params->svd_paramslist[range_loop_counter].n_range_end;
                q_inc = params->svd_paramslist[range_loop_counter].n_range_step_size;
                params->datatype = params->svd_paramslist[range_loop_counter].data_types;
                params->datatype_char = params->svd_paramslist[range_loop_counter].data_types_char;
                n_repeats = params->svd_paramslist[range_loop_counter].num_repeats;
                n_datatypes = params->svd_paramslist[range_loop_counter].num_data_types;
                break;

            case AUX:
                p_first = params->aux_paramslist[range_loop_counter].m_range_start;
                p_max = params->aux_paramslist[range_loop_counter].m_range_end;
                p_inc = params->aux_paramslist[range_loop_counter].m_range_step_size;
                q_first = params->aux_paramslist[range_loop_counter].n_range_start;
                q_max = params->aux_paramslist[range_loop_counter].n_range_end;
                q_inc = params->aux_paramslist[range_loop_counter].n_range_step_size;
                params->datatype = params->aux_paramslist[range_loop_counter].data_types;
                params->datatype_char = params->aux_paramslist[range_loop_counter].data_types_char;
                n_repeats = params->aux_paramslist[range_loop_counter].num_repeats;
                n_datatypes = params->aux_paramslist[range_loop_counter].num_data_types;
                break;

            default:
                return;
        }

        /* Loop over the requested datatypes. */
        for(dt = 0; dt < n_datatypes; ++dt)
        {
            datatype = params->datatype[dt];
            datatype_char = params->datatype_char[dt];
            /* Skip complex and double complex tests of not supported APIs */
            if(!FLA_SKIP_TEST(datatype_char, func_str))
            {
                /* Loop over the requested problem sizes */
                for(p_cur = p_first, q_cur = q_first; (p_cur <= p_max && q_cur <= q_max);
                    p_cur += p_inc, q_cur += q_inc)
                {
                    params->n_repeats = n_repeats;
                    if(n_threads > 1)
                    {
#pragma omp parallel num_threads(n_threads)
#pragma omp for
                        for(ith = 0; ith < n_threads; ith++)
                        {
                            f_exp(func_str, params, datatype, p_cur, q_cur, range_loop_counter,
                                  n_repeats, einfo);
                        }
                    }
                    else
                    {
                        f_exp(func_str, params, datatype, p_cur, q_cur, range_loop_counter,
                              n_repeats, einfo);
                    }
                }
            }
        }

        fla_test_output_info("\n");
    }

    free(perf);
    free(time);
}

/* This function is used to clock the execution time of APIs*/
double fla_test_clock()
{
#ifdef _WIN32
    LARGE_INTEGER clock_freq = {0};
    LARGE_INTEGER clock_val;

    QueryPerformanceFrequency(&clock_freq);
    QueryPerformanceCounter(&clock_val);

    return ((double)clock_val.QuadPart / (double)clock_freq.QuadPart);
#else
    double the_time, norm_sec;
    struct timespec ts;

    clock_gettime(CLOCK_MONOTONIC, &ts);

    if(ref_time_sec == 0.0)
        ref_time_sec = (double)ts.tv_sec;

    norm_sec = (double)ts.tv_sec - ref_time_sec;

    the_time = norm_sec + ts.tv_nsec * 1.0e-9;

    return the_time;
#endif
}

/*
 * This function checks if the argument "--bench=<k>" is present in the command line.
 * If the argument is present, it sets the benchmark_mode flag and bench_duration in
 * params.
 *
 * In the benchmark mode, api execution would be run atleast for <k> seconds.
 * It also shows more detailed statistics by default.
 *
 * @return
 *      0 if the argument is not found,
 *      1 if the argument is found and valid,
 *     -1 if the argument is found but invalid value is provided.
 * */
integer fla_parse_bench_arg(integer argc, char **argv, test_params_t *params)
{
    params->benchmark_mode = 0;
    const char *arg_str = "--bench=";
    integer arg_len = strlen(arg_str);

    /* Check if the argument is present in the command line. */
    for(integer i = 1; i < argc; ++i)
    {
        if(strncmp(argv[i], arg_str, arg_len) == 0)
        {
            errno = 0; /* Reset errno before strtof */
            params->bench_duration = strtof(argv[i] + arg_len, NULL);
            /* Checking any error during parsing */
            /* Checking if bench value is not negative */
            if(errno != 0 || params->bench_duration <= 0)
            {
                /* Invalid bench value */
                printf("\nError: Invalid bench argument: %s\n", argv[i]);
                printf("       Please provide a valid duration > 0.\n");
                return -1;
            }
            /* Bench value is valid */
            params->benchmark_mode = 1;
            return 1;
        }
    }
    return 0;
}

/*
 * This function checks if the argument "--warmup=<k>" is present in the command line.
 * If the argument is present, it sets the warmup_repeats in params.
 * The value can be:
 *      0: Disable Warmup
 *      Any decimal k in range (0, 1) exclusive: ceil(k * repeat) number of invocations as warmup
 *      Any integer k >= 1: k invocations as warmup
 * @return
 *      0 if the argument is not found,
 *      1 if the argument is found and valid,
 *     -1 if the argument is found but invalid value is provided.
 */
integer fla_parse_warmup_arg(integer argc, char **argv, test_params_t *params)
{
    const char *arg_str = "--warmup=";
    integer arg_len = strlen(arg_str);

    /* Check if the argument is present in the command line. */
    for(integer i = 1; i < argc; ++i)
    {
        if(strncmp(argv[i], arg_str, arg_len) == 0)
        {
            errno = 0; /* Reset errno before strtof */
            params->warmup_repeats = strtof(argv[i] + arg_len, NULL);
            /* Checking any error during parsing */
            /* Checking if warmup value is not negative */
            /* Checking if warmup value is > 1, then it is an integer */
            if(errno != 0 || params->warmup_repeats < 0
               || (params->warmup_repeats > 1
                   && (params->warmup_repeats - (integer)params->warmup_repeats) != 0))
            {
                /* Invalid warmup value */
                printf("\nError: Invalid warmup argument: %s\n", argv[i]);
                printf("       Please provide a valid number.\n");
                printf("       0: Disable Warmup\n");
                printf("       Any decimal k in range (0, 1) exclusive: "
                       "ceil(k * repeat) number of invocations as warmup\n");
                printf("       Any integer k >= 1: k invocations as warmup\n");
                return -1;
            }
            /* Warmup value is valid */
            return 1;
        }
    }
    return 0;
}

/*
 * This function checks if the argument "--stats=<stats_list>" is present in the command line.
 * stats_list is a comma separated list of statistics to be printed.
 * The list can contain the following values:
 *      min: Minimum time
 *      max: Maximum time
 *      avg: Average time
 *      p1-p99: Percentiles
 * The function sets the stats_out array in params with the corresponding values.
 * @return
 *      0 if the argument is not found,
 *      1 if the argument is found and valid,
 *     -1 if the argument is found but invalid value is provided.
 */
integer fla_parse_stats_arg(integer argc, char **argv, test_params_t *params)
{
    params->num_stats = 0;
    const char *arg_str = "--stats=";
    integer arg_len = strlen(arg_str);
    char *stats_val = NULL;
    integer num_stats_found = 0;
    integer percentile_val;
    char *endptr = NULL;

    /* Check if the argument is present in the command line. */
    for(integer i = 1; i < argc; ++i)
    {
        if(strncmp(argv[i], arg_str, arg_len) == 0)
        {
            stats_val = argv[i] + arg_len;
            /* Set values to lowercase */
            for(integer j = 0; j < strlen(stats_val); j++)
            {
                stats_val[j] = tolower(stats_val[j]);
            }
            /* Split by comma */
            char *token = stats_val;
            char *token_next = NULL;
            integer token_len = 0;
            integer found_flag;
            while(token != NULL)
            {
                token_next = strchr(token, ',');

                /* This is the last token */
                if(token_next == NULL)
                {
                    token_len = strlen(token);
                }
                else
                {
                    token_len = token_next - token;
                }

                if(token_len <= 0)
                {
                    num_stats_found = 0;
                    break;
                }

                if(num_stats_found >= MAX_NUM_STATS)
                {
                    printf("\nError: Too many stats specified. Max allowed is %d\n", MAX_NUM_STATS);
                    return -1;
                }

                found_flag = 0;
                for(integer j = 0; j < FLA_NUM_STATS; j++)
                {
                    integer stat_str_len = strlen(AVAILABLE_STATS[j].arg_str);
                    if(strncmp(token, AVAILABLE_STATS[j].arg_str, stat_str_len) == 0)
                    {
                        /* If percentile then need special input */
                        if(AVAILABLE_STATS[j].stat_type == FLA_PERCENTILE_STAT)
                        {
                            if((token_len - stat_str_len) >= 1 && (token_len - stat_str_len) <= 2)
                            {
                                percentile_val = (integer)strtol(token + stat_str_len, &endptr, 10);
                                if((*endptr == '\0' || *endptr == ',') && percentile_val >= 1
                                   && percentile_val <= 99)
                                {
                                    params->stats_out[num_stats_found++]
                                        = (fla_stat_t){.stat_type = AVAILABLE_STATS[j].stat_type,
                                                       .percentile_num = percentile_val};
                                    found_flag = 1;
                                }
                            }
                        }
                        else if(token_len == stat_str_len)
                        {
                            params->stats_out[num_stats_found++].stat_type
                                = AVAILABLE_STATS[j].stat_type;
                            found_flag = 1;
                        }
                    }
                }
                if(found_flag == 0)
                {
                    num_stats_found = 0;
                    break;
                }

                token = token_next != NULL ? token_next + 1 : NULL;
            }

            if(num_stats_found == 0)
            {
                printf("\nError: Invalid stats argument: %s\n", argv[i]);
                printf("       Please provide a valid stats argument.\n");
                printf("       Upto %d stats can be specified separated by comma.\n",
                       MAX_NUM_STATS);
                printf("       Available stats are:");
                for(integer j = 0; j < FLA_NUM_STATS - 1; j++)
                {
                    if(AVAILABLE_STATS[j].stat_type == FLA_PERCENTILE_STAT)
                    {
                        printf(" %s<1-99>,", AVAILABLE_STATS[j].arg_str);
                    }
                    else
                    {
                        printf(" %s,", AVAILABLE_STATS[j].arg_str);
                    }
                }
                printf(" %s\n", AVAILABLE_STATS[FLA_NUM_STATS - 1].arg_str);
                printf("       Example: --stats=min,avg,p95,var\n");
                return -1;
            }

            params->num_stats = num_stats_found;
            return 1;
        }
    }
    return 0;
}

/*
 * This function checks if the argument "--print-header" is present in the command line.
 * --print-header is used to print the header for the test output in cli mode.
 *
 * * @return
 *      0 if the argument is not found,
 *      1 if the argument is found .
 **/
integer fla_parse_print_header_arg(integer argc, char **argv, test_params_t *params)
{
    const char *arg_str = "--print-header";
    integer arg_len = strlen(arg_str);

    /* Check if the argument is present in the command line. */
    for(integer i = 1; i < argc; ++i)
    {
        if(strncmp(argv[i], arg_str, arg_len) == 0)
        {
            params->cli_print_header = 1;
            return 1;
        }
    }
    return 0;
}

/*
 * This function parses the time unit argument "--time-unit=<unit>" from the command line.
 * The unit can be one of the following: s, ms, us, ns, ps, auto.
 *
 * * @return
 *      0 if the argument is not found,
 *      1 if the argument is found and valid,
 *     -1 if the argument is found but invalid value is provided.
 */
integer fla_parse_time_unit_arg(integer argc, char **argv, test_params_t *params)
{
    const char *arg_str = "--time-unit=";
    integer arg_len = strlen(arg_str);
    char *unit_arg = NULL;

    /* Check if the argument is present in the command line. */
    for(integer i = 1; i < argc; ++i)
    {
        if(strncmp(argv[i], arg_str, arg_len) == 0)
        {
            unit_arg = argv[i] + arg_len;
            /* Set values to lowercase */
            for(integer j = 0; j < strlen(unit_arg); j++)
            {
                unit_arg[j] = tolower(unit_arg[j]);
            }
            if(strncmp(unit_arg, "s", 1) == 0)
            {
                params->time_unit = FLA_TIME_UNIT_SEC;
            }
            else if(strncmp(unit_arg, "ms", 2) == 0)
            {
                params->time_unit = FLA_TIME_UNIT_MILLISEC;
            }
            else if(strncmp(unit_arg, "us", 2) == 0)
            {
                params->time_unit = FLA_TIME_UNIT_MICROSEC;
            }
            else if(strncmp(unit_arg, "ns", 2) == 0)
            {
                params->time_unit = FLA_TIME_UNIT_NANOSEC;
            }
            else if(strncmp(unit_arg, "ps", 2) == 0)
            {
                params->time_unit = FLA_TIME_UNIT_PICOSEC;
            }
            else if(strncmp(unit_arg, "auto", 4) == 0)
            {
                params->time_unit = FLA_TIME_UNIT_AUTO;
            }
            else
            {
                printf("\nError: Invalid time unit argument: %s\n", argv[i]);
                printf("       Please provide a valid time unit argument.\n");
                printf("       Available time units are: s, ms, us, ns, ps, auto\n");
                return -1;
            }
            return 1;
        }
    }
    return 0;
}

/*
 * This function checks if the argument "--help" is present in the command line.
 * If the argument is present, it prints the help message and returns -1.
 *
 * * @return
 *     0 if the argument is not found,
 *    -1 if the argument is found and help message is printed.
 **/
integer fla_check_help_arg(integer argc, char **argv)
{
    for(integer i = 1; i < argc; ++i)
    {
        if(strcmp(argv[i], "--help") == 0)
        {
            fla_print_help(argv[0]);
            return -1;
        }
    }
    return 0;
}

/*
 * This function checks if the argument "--filter-outliers[=<multiplier>]" is present in the command
 * line. If the argument is present, it will filter runtimes greater than
 *
 *   multiplier * stddev + mean
 *
 * If multiplier is not provided, it will use the default value of FLA_OUTLIERS_MULTIPLIER_DEFAULT.
 *
 * @return
 *      0 if the argument is not found,
 *      1 if the argument is found and valid,
 *     -1 if the argument is found but invalid value is provided.
 */
integer fla_parse_filter_outliers_args(integer argc, char **argv, test_params_t *params)
{
    const char *arg_str = "--filter-outliers";
    integer arg_len = strlen(arg_str);

    /* Check if the argument is present in the command line. */
    for(integer i = 1; i < argc; ++i)
    {
        if(strncmp(argv[i], arg_str, arg_len) == 0)
        {
            /* If '=' is present then filter parameter is provided */
            if(argv[i][arg_len] == '=')
            {
                errno = 0; /* Reset errno before strtof */
                params->outlier_multiplier = strtof(argv[i] + arg_len + 1, NULL);
                /* Checking any error during parsing */
                /* Checking if filter value is not negative */
                if(errno != 0 || params->outlier_multiplier <= 0)
                {
                    /* Invalid filter value */
                    printf("\nError: Invalid filter-outliers argument: %s\n", argv[i]);
                    printf("       Please provide a valid multiplier > 0.\n");
                    return -1;
                }
            }
            else
            {
                /* If '=' is not present then default muliplier is used */
                params->outlier_multiplier = FLA_OUTLIERS_MULTIPLIER_DEFAULT;
            }
            return 1;
        }
    }
    return 0;
}

/* This function checks if the argument "--dump-runtimes=<file_name>" is present in the command
 * line. If the argument is present, it sets the dump_runtimes_file_name in params.
 *
 * @return
 *      0 if the argument is not found,
 *      1 if the argument is found and valid,
 *     -1 if the argument is found but invalid value is provided.
 */
integer fla_parse_dump_runtimes_arg(integer argc, char **argv, test_params_t *params)
{
    const char *arg_str = "--dump-runtimes=";
    integer arg_len = strlen(arg_str);

    /* Check if the argument is present in the command line. */
    for(integer i = 1; i < argc; ++i)
    {
        if(strncmp(argv[i], arg_str, arg_len) == 0)
        {
            if(argv[i][arg_len] == '\0')
            {
                printf("\nError: Invalid dump-runtimes argument: %s\n", argv[i]);
                printf("       Please provide a valid file name.\n");
                return -1;
            }
            else
            {
                params->dump_runtimes_file_name = argv[i] + arg_len;
                return 1;
            }
        }
    }
    return 0;
}

#define FLA_ARGS_PARSE_RESULT_HANDLER                 \
    if(parse_status < 0)                              \
    {                                                 \
        return FALSE; /* Invalid arg while parsing */ \
    }                                                 \
    n_args_found += parse_status;

/*
 * This function parses the command line arguments and sets the
 * corresponding values in the params structure.
 * It also updates the argc values to ignore parsed arguments in further
 * processing. These arguments should be passed in the last.
 *
 * @return TRUE if parsing is successful, FALSE if any error occured
 * */
bool fla_parse_cmdline_args(integer *argc, char **argv, test_params_t *params)
{
    integer n_args_found = 0;
    integer parse_status;

    /* Check if the help argument is present */
    parse_status = fla_check_help_arg(*argc, argv);
    FLA_ARGS_PARSE_RESULT_HANDLER;

    parse_status = fla_check_interface(*argc, argv, params);
    FLA_ARGS_PARSE_RESULT_HANDLER;

    parse_status = fla_parse_bench_arg(*argc, argv, params);
    FLA_ARGS_PARSE_RESULT_HANDLER;

    parse_status = fla_parse_warmup_arg(*argc, argv, params);
    FLA_ARGS_PARSE_RESULT_HANDLER;

    /* If warmup is not provided and benchmark mode then set
    default warmup */
    if(parse_status == 0 && params->benchmark_mode)
    {
        params->warmup_repeats = FLA_BENCH_DEFAULT_WARMUP;
    }

    parse_status = fla_parse_stats_arg(*argc, argv, params);
    FLA_ARGS_PARSE_RESULT_HANDLER;
    /* If stats arg not found then set default stats */
    if(params->num_stats == 0)
    {
        if(params->benchmark_mode)
        {
            /* In benchmark mode show min, avg and p95 */
            params->stats_out[0] = (fla_stat_t){.stat_type = FLA_MIN_STAT};
            params->stats_out[1] = (fla_stat_t){.stat_type = FLA_AVG_STAT};
            params->stats_out[2]
                = (fla_stat_t){.stat_type = FLA_PERCENTILE_STAT, .percentile_num = 95};
            params->num_stats = 3;
        }
        else
        {
            /* In normal mode show only min */
            params->stats_out[0] = (fla_stat_t){.stat_type = FLA_MIN_STAT};
            params->num_stats = 1;
        }
    }

    parse_status = fla_parse_print_header_arg(*argc, argv, params);
    FLA_ARGS_PARSE_RESULT_HANDLER;

    parse_status = fla_parse_time_unit_arg(*argc, argv, params);
    FLA_ARGS_PARSE_RESULT_HANDLER;

    parse_status = fla_parse_filter_outliers_args(*argc, argv, params);
    FLA_ARGS_PARSE_RESULT_HANDLER;

    parse_status = fla_parse_dump_runtimes_arg(*argc, argv, params);
    FLA_ARGS_PARSE_RESULT_HANDLER;

    *argc -= n_args_found;

    return TRUE;
}

/**
 * This function checks if the runtimes array is needed.
 *
 * @return 1 if runtimes array is needed, 0 otherwise.
 */
integer fla_need_runtimes_array(test_params_t *params)
{
    if(params->dump_runtimes_file_name || params->outlier_multiplier)
    {
        return 1;
    }
    for(integer i = 0; i < params->num_stats; i++)
    {
        if(params->stats_out[i].stat_type == FLA_PERCENTILE_STAT
           || params->stats_out[i].stat_type == FLA_VARIANCE_STAT
           || params->stats_out[i].stat_type == FLA_STDDEV_STAT)
        {
            return 1;
        }
    }
    return 0;
}

/* Utility functions to get the statistics values
 * If percentile is needed then make sure filtered_runarr is sorted.
 */
double fla_stat_get_val(test_params_t *params, fla_stat_t *stat_desc)
{
    test_runtime_ctx_t *ctx = &params->runtime_ctx;

    switch(stat_desc->stat_type)
    {
        case FLA_MIN_STAT:
        {
            double min_time = ctx->min_time;
            if(params->outlier_multiplier != 0.0 && ctx->run_times_arr != NULL)
            {
                get_min_from_array(DOUBLE, ctx->run_times_arr, &min_time,
                                   ctx->filtered_run_times_size);
            }
            return min_time;
        }
        case FLA_MAX_STAT:
        {
            double max_time = ctx->max_time;
            if(params->outlier_multiplier != 0.0 && ctx->run_times_arr != NULL)
            {
                get_max_from_array(DOUBLE, ctx->run_times_arr, &max_time,
                                   ctx->filtered_run_times_size);
            }
            return max_time;
        }
        case FLA_AVG_STAT:
        {
            double avg = ctx->total_time / ctx->run_times_counter;
            if(params->outlier_multiplier != 0.0 && ctx->run_times_arr != NULL)
            {
                get_avg_of_array(DOUBLE, ctx->run_times_arr, &avg, ctx->filtered_run_times_size);
            }
            return avg;
        }

        case FLA_PERCENTILE_STAT:
        {
            size_t p_index
                = (size_t)ceil((ctx->filtered_run_times_size * stat_desc->percentile_num) / 100.0);
            p_index = fla_min(p_index, ctx->filtered_run_times_size - 1);
            return ctx->run_times_arr != NULL ? ctx->run_times_arr[p_index] : -1.0;
        }
        case FLA_VARIANCE_STAT:
        {
            double variance = 0.0;
            if(ctx->run_times_arr != NULL)
            {
                get_variance_of_array(DOUBLE, ctx->run_times_arr, &variance,
                                      ctx->filtered_run_times_size);
            }
            return variance;
        }
        case FLA_STDDEV_STAT:
        {
            double stddev = 0.0;
            if(ctx->run_times_arr != NULL)
            {
                get_stddev_of_array(DOUBLE, ctx->run_times_arr, &stddev,
                                    ctx->filtered_run_times_size);
            }
            return stddev;
        }
        default:
        {
            return -1.0;
        }
    }
    /* If none of the above cases are met, return -1.0 */
    return -1.0;
}

/* This function initializes the runtime context */
void fla_test_runtime_ctx_init(test_params_t *params)
{
    params->runtime_ctx.run_times_counter = 0;
    params->runtime_ctx.run_times_arr = NULL;
    params->runtime_ctx.run_times_arr_size = 0;
    params->runtime_ctx.warmup_counter = 0;
    params->runtime_ctx.min_time = DBL_MAX;
    params->runtime_ctx.total_time = 0.0;
    params->runtime_ctx.max_time = 0.0;
    params->runtime_ctx.filtered_run_times_size = 0;
    params->runtime_ctx.need_runtimes_array = fla_need_runtimes_array(params);
}

/* Resets the runtime context. Does not free the run_times_arr
 * as this array can be used for multiple tests.
 */
void fla_test_runtime_ctx_reset(test_params_t *params)
{
    params->runtime_ctx.run_times_counter = 0;
    params->runtime_ctx.warmup_counter = 0;
    params->runtime_ctx.min_time = DBL_MAX;
    params->runtime_ctx.max_time = 0.0;
    params->runtime_ctx.total_time = 0.0;
    params->runtime_ctx.filtered_run_times_size = 0;
}

/* Frees memory allocated for the context.
 * This function should be called when all tests are finished
 */
void fla_test_runtime_ctx_free(test_params_t *params)
{
    if(params->runtime_ctx.run_times_arr != NULL && params->runtime_ctx.run_times_arr_size > 0)
    {
        free(params->runtime_ctx.run_times_arr);
        params->runtime_ctx.run_times_arr = NULL;
        params->runtime_ctx.run_times_arr_size = 0;
    }
}