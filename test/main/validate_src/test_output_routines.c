/******************************************************************************
 * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

#include "test_common.h"
#include "test_lapack.h"
#include "test_output_routines.h"

#define MAX_LINE_LENGTH 256

void fla_test_output_info(char *message, ...)
{
    FILE *output_stream = stdout;
    va_list args;

    /* Initialize variable argument environment. */
    va_start(args, message);

    /* Parse the received message and print its components. */
    vfprintf(output_stream, message, args);

    /* Shutdown variable argument environment and clean up stack. */
    va_end(args);

    /* Flush the output stream. */
    fflush(output_stream);
}

void fla_test_output_error(char *message, ...)
{
    FILE *output_stream = stderr;
    va_list args;

    fprintf(output_stream, "%s: *** error ***: ", fla_test_binary_name);

    // Initialize variable argument environment.
    va_start(args, message);

    // Parse the received message and print its components.
    vfprintf(output_stream, message, args);

    // Shutdown variable argument environment and clean up stack.
    va_end(args);

    // Flush the output stream.
    fflush(output_stream);

    // Exit.
    exit(1);
}

#define FLA_EXIT_IF_STR_LEN_EXCEEDS(STRLEN, MAXLEN) \
    if(STRLEN >= MAXLEN)                            \
    {                                               \
        return;                                     \
    }

void fla_test_print_header(void *test_params)
{
    char header[MAX_LINE_LENGTH];
    char underline[MAX_LINE_LENGTH];

    integer header_i = 0;
    integer underline_i = 0;
    integer write_len;
    int space_len;

    test_params_t *params = (test_params_t *)test_params;

    header_i += snprintf(header, MAX_LINE_LENGTH, "%2sAPI%13s DATA_TYPE%6s SIZE%9s FLOPS%7s ", "",
                         "", "", "", "");

    FLA_EXIT_IF_STR_LEN_EXCEEDS(header_i, MAX_LINE_LENGTH);

    /* If benchmark mode than add the number of repeats that actually happened */
    if(params->benchmark_mode)
    {
        header_i += snprintf(header + header_i, MAX_LINE_LENGTH - header_i, "REPEATS%7s ", "");
        FLA_EXIT_IF_STR_LEN_EXCEEDS(header_i, MAX_LINE_LENGTH);
    }

    /* In case of manual time unit, add extra space to the time column
       becasue scientific notation is used for time values */
    char *time_extra_space = params->time_unit == FLA_TIME_UNIT_AUTO ? "" : "       ";

    /* If only min stat is required in normal mode, then just print TIME */
    if(!params->benchmark_mode && params->num_stats == 1
       && params->stats_out[0].stat_type == FLA_MIN_STAT)
    {
        header_i += snprintf(header + header_i, MAX_LINE_LENGTH - header_i, "  TIME%9s %s",
                              "", time_extra_space);
        FLA_EXIT_IF_STR_LEN_EXCEEDS(header_i, MAX_LINE_LENGTH);
    }
    else
    {
        for(integer i = 0; i < params->num_stats; ++i)
        {
            if(params->stats_out[i].stat_type == FLA_PERCENTILE_STAT)
            {
                write_len
                    = snprintf(header + header_i, MAX_LINE_LENGTH - header_i, "TIME(%s%" FT_IS ")%s",
                               AVAILABLE_STATS[params->stats_out[i].stat_type].out_str,
                               params->stats_out[i].percentile_num, time_extra_space);
            }
            else
            {
                write_len = snprintf(header + header_i, MAX_LINE_LENGTH - header_i, "TIME(%s)%s",
                                     AVAILABLE_STATS[params->stats_out[i].stat_type].out_str,
                                     time_extra_space);
            }
            header_i += write_len;
            FLA_EXIT_IF_STR_LEN_EXCEEDS(header_i, MAX_LINE_LENGTH);
            space_len = fla_max(14 - write_len, 0);
            header_i
                += snprintf(header + header_i, MAX_LINE_LENGTH - header_i, "%*s", space_len, "");
            FLA_EXIT_IF_STR_LEN_EXCEEDS(header_i, MAX_LINE_LENGTH);
        }
    }

    header_i += snprintf(header + header_i, MAX_LINE_LENGTH - header_i, "ERROR%9s STATUS", "");
    FLA_EXIT_IF_STR_LEN_EXCEEDS(header_i, MAX_LINE_LENGTH);

    /* Filling up underline */
    integer underline_len = fla_min(header_i + 1, MAX_LINE_LENGTH - 1);
    for(underline_i = 0; underline_i < underline_len; ++underline_i)
    {

        if((underline_i < (header_i - 1) && header[underline_i + 1] != ' ')
           || (underline_i > 0 && header[underline_i - 1] != ' '))
        {
            underline[underline_i] = '=';
        }
        else
        {
            underline[underline_i] = ' ';
        }
    }
    underline[underline_i] = '\0';

    /* Print the header */
    fla_test_output_info("%s\n", header);
    fla_test_output_info("%s\n", underline);
}

void fla_print_help(char *exe_name)
{
    printf("Usage:\n\n");
    printf("With config file: \n%s [--config-dir=<config_dir_path>] [test_suite_options]\n",
           exe_name);
    printf("\n  <config_dir_path> is the path to the config directory. The following \n"
           "  configurations are provided by default based on input problem sizes: micro, short, "
           "medium, long\n");
    printf("  If --config-dir is not specified, the short config directory is used by default.\n");
    printf(
        "\nWith command line args: \n%s <api_name> <api_args> <repeats> [input_validation_options] "
        "[test_suite_options]\n",
        exe_name);
    printf("\n  <api_name> is the name of the API to be tested and <api_args> are API-specific "
           "test parameters\n");
    printf("  <repeats> is the number of times the test should be run.\n");

    printf("\nInput validation options (Optional):\n");
    printf("Only one of the following options can be used at a time.\n");

    printf("\n  --imatrix=<type>: Use special extreme input values. <type> can be one of the "
           "following:\n");
    printf("       O: Perform overflow test\n");
    printf("       U: Perform underflow test\n");
    printf("       N: Initialize the matrix with NAN values in all locations\n");
    printf("       I: Initialize the matrix with INFINITY values in all locations\n");
    printf("       A: Initialize the matrix with NAN values in a few random locations\n");
    printf("       F: Initialize the matrix with INFINITY values in a few random locations\n");

    printf("\n  --einfo=<err_code>: Validate that the API returns the expected error code.\n");

    printf("\n  <file_path>: Path to the input file containing the input data. Input matrices are "
           "populated from this file\n");

    printf("\nTest suite options (Optional):\n");

    printf("\n  --interface=<interface>: Use the specified interface for the test.");
    printf("\n  If both --interface and --lapacke are specified, --interface will be ignored.\n");
    printf("  <interface> can be one of the following:\n");
    printf("       lapack: Use LAPACK interface\n");
    printf("       cpp: Use C++ interface\n");
    printf("       lapacke_row: Use LAPACKE interface with row-ordered matrices\n");
    printf("       lapacke_column: Use LAPACKE interface with column-ordered matrices\n");

    printf("\n  --bench=<k>: Run the test in benchmark mode. <k> is the duration in seconds for "
           "which the "
           "test should be run.\n");
    printf("  The actual number of repeats is calculated as max(repeats, total_test_time)\n");
    printf("  In benchmark mode, by default %f%% of invocations are run as warmup.\n",
           FLA_BENCH_DEFAULT_WARMUP);

    printf("\n  --warmup=<k>: Number of warmup invocations. <k> can be:\n");
    printf("       0: Disable warmup\n");
    printf("       Any decimal k in range (0, 1) exclusive: ceil(k * repeat) number of invocations "
           "as warmup\n");
    printf("       Any integer k >= 1: k invocations as warmup\n");

    printf("\n  --stats=<stats_list>: Comma-separated list of statistics to be printed.\n");
    printf("       Up to %d stats can be specified.\n", MAX_NUM_STATS);
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
    printf("       If not specified, default stats are min, avg, and p95 in benchmark mode.\n");
    printf("       In normal mode, only min is shown.\n");

    printf("\n  --print-header: Print the header for the test output in CLI mode.\n");
    printf("       By default, the header is not printed in CLI mode.\n");

    printf("\n  --time-unit=<unit>: Specify the time unit for the test output.\n");
    printf("       <unit> can be one of the following: s, ms, us, ns, ps, auto.\n");
    printf("       If not specified, the default unit is auto.\n");
    printf("       In auto mode, the time unit is automatically selected based on the time taken "
           "for the test.\n");

    printf(
        "\n  --filter-outliers[=<multiplier>]: Filter out the outliers from the test results.\n");
    printf("       If <multiplier> is not specified, the default value is %lf.\n",
           FLA_OUTLIERS_MULTIPLIER_DEFAULT);
    printf("       Filters out the values greater than <multiplier> * stddev + mean.\n");

    printf("\n  --dump-runtimes=<file_path>: Dump the runtimes to the specified file.\n");
    printf("       This option is only valid in CLI mode.\n");

    printf("\n  --help: Show this help message.\n");
}

void fla_test_print_status(char *func_str, char datatype_char, integer sqr_inp, integer p_cur,
                           integer q_cur, double residual, double thresh, double perf,
                           void *test_params)
{
    integer datatype, i;
    char *pass_str;
    double scale_factor = 1.0;
    char t_scale[3] = "";
    double t_scale_factor = 1.0;

    test_params_t *params = (test_params_t *)test_params;
    size_t eff_repeats = params->runtime_ctx.run_times_counter;

    double stats_vals[params->num_stats];
    char stats_scales[params->num_stats][3];
    datatype = get_datatype(datatype_char);

    FLA_MAP_API_NAME(datatype, func_str);

    pass_str = fla_test_get_string_for_result(residual, datatype, thresh);

    /* Populate stats val */
    if(fla_populate_stat_vals(stats_vals, params) != 0)
    {
        fla_test_output_error("Error: Failed to populate stats values\n");
        return;
    }

    /* Decide time unit based on the maximum value */
    double max_stats_val;
    get_max_from_array(DOUBLE, stats_vals, &max_stats_val, params->num_stats);

    fla_test_get_time_unit(t_scale, max_stats_val, &scale_factor, params->time_unit);

    /* Convert the time to the required unit */
    for(i = 0; i < params->num_stats; ++i)
    {
        if(params->stats_out[i].stat_type != FLA_VARIANCE_STAT
           && params->stats_out[i].stat_type != FLA_STDDEV_STAT)
        {
            stats_vals[i] *= scale_factor;
            strcpy(stats_scales[i], t_scale);
        }
        else
        {
            /* Since variance and stddev can be very small, scale them separately */
            fla_test_get_time_unit(stats_scales[i], stats_vals[i], &t_scale_factor,
                                   FLA_TIME_UNIT_AUTO);
            stats_vals[i] *= t_scale_factor;
        }
    }

    if(sqr_inp == SQUARE_INPUT)
    {
        q_cur = p_cur;
    }

    fla_test_output_info(" %-20s  %c  %10" FT_IS " x %-9" FT_IS " %-10.2lf ", func_str,
                         datatype_char, p_cur, q_cur, perf);

    if(params->benchmark_mode)
    {
        fla_test_output_info(" %9zu %4s", eff_repeats, "");
    }

    /* In case of manual time unit, use scientific notation for time values */
    if(!params->benchmark_mode && params->num_stats == 1
       && params->stats_out[0].stat_type == FLA_MIN_STAT)
    {
        fla_test_output_info(params->time_unit == FLA_TIME_UNIT_AUTO ? " %6.2lf %-7s " : " %6.2e %-7s ", stats_vals[0], stats_scales[0]);
    }
    else
    {
        for(i = 0; i < params->num_stats; ++i)
        {
            fla_test_output_info(params->time_unit == FLA_TIME_UNIT_AUTO ? " %6.2lf %-5s " : " %6.2e %-5s ", stats_vals[i], stats_scales[i]);
        }
    }

    fla_test_output_info(" %-7.2le ", residual);

    if(pass_str[0] == 'P' || pass_str[0] == 'F')
    {
        fla_test_output_info("  %10s\n", pass_str);
    }
    else
    {
        fla_test_output_info("  %12s\n", pass_str);
    }

    g_total_tests++;
    g_tests_passed[datatype - FLOAT] += (pass_str[0] == 'P');
    g_tests_failed[datatype - FLOAT] += (pass_str[0] == 'F');
    g_tests_incomplete[datatype - FLOAT] += (pass_str[0] == 'I');
    g_total_failed_tests += (pass_str[0] == 'F');
    g_total_incomplete_tests += (pass_str[0] == 'I');
}

/* This function sets the unit of time
 * and provides the scale factor to convert the time
 * from seconds to the required unit.
 */
void fla_test_get_time_unit(char *scale, double measured_time, double *scale_factor,
                            integer time_unit)
{
    switch(time_unit)
    {
        case FLA_TIME_UNIT_AUTO:
        {
            /* Automatically determine appropriate time unit */
            if(measured_time >= 1)
            {
                strcpy(scale, "s");
                *scale_factor = 1;
                return;
            }

            if((measured_time < 1) && (measured_time >= 0.001))
            {
                strcpy(scale, "ms");
                *scale_factor = 1000;
            }
            else if((measured_time < 0.001) && (measured_time >= 0.000001))
            {
                strcpy(scale, "us");
                *scale_factor = 1000000;
            }
            else if((measured_time < 0.000001) && (measured_time >= 0.000000001))
            {
                strcpy(scale, "ns");
                *scale_factor = 1000000000;
            }
            else if((measured_time < 0.000000001))
            {
                strcpy(scale, "ps");
                *scale_factor = 1000000000000;
            }
            break;
        }
        /* Fixed time unit required */
        case FLA_TIME_UNIT_SEC:
        {
            strcpy(scale, "s");
            *scale_factor = 1;
            break;
        }
        case FLA_TIME_UNIT_MILLISEC:
        {
            strcpy(scale, "ms");
            *scale_factor = 1000;
            break;
        }
        case FLA_TIME_UNIT_MICROSEC:
        {
            strcpy(scale, "us");
            *scale_factor = 1000000;
            break;
        }
        case FLA_TIME_UNIT_NANOSEC:
        {
            strcpy(scale, "ns");
            *scale_factor = 1000000000;
            break;
        }
        case FLA_TIME_UNIT_PICOSEC:
        {
            strcpy(scale, "ps");
            *scale_factor *= 1000000000000;
            break;
        }
        default:
        {
            fla_test_output_error("Error: Invalid time unit %d\n", time_unit);
            break;
        }
    }
}

/*
 * Populate the stat_vals array with the required statistics
 * returns 0 on success, -1 on failure
 */
integer fla_populate_stat_vals(double *stat_vals, void *test_params)
{
    if(stat_vals == NULL)
    {
        fla_test_output_error("Error: stat_vals is NULL\n");
        return -1;
    }

    test_params_t *params = (test_params_t *)test_params;
    test_runtime_ctx_t *ctx = &params->runtime_ctx;

    integer i;
    integer is_sort_needed = 0;

    fla_filter_runtimes(params);

    if(ctx->need_runtimes_array && ctx->filtered_run_times_size == 0)
    {
        for(i = 0; i < params->num_stats; ++i)
        {
            stat_vals[i] = 0.0;
        }
        fla_test_output_info("Warning: No valid runtimes found for the test after filtering.\n"
                             "Setting all stats to 0.0\n");
        return 0;
    }

    /* IF Percentile specified then sort is needed */
    for(i = 0; i < params->num_stats; ++i)
    {
        if(params->stats_out[i].stat_type == FLA_PERCENTILE_STAT)
        {
            is_sort_needed = 1;
            break;
        }
    }

    if(is_sort_needed && ctx->filtered_run_times_size != 0)
    {
        /* First sort the array */
        qsort_realtype_vector(DOUBLE, "A", ctx->filtered_run_times_size, ctx->run_times_arr);
    }

    for(i = 0; i < params->num_stats; ++i)
    {
        stat_vals[i] = fla_stat_get_val(params, &params->stats_out[i]);
    }

    return 0;
}

void fla_test_print_summary()
{
    fla_test_output_info("\n\nResults Summary:\n\n");
    fla_test_output_info("%2sDATATYPE%13s No. of Tests%6s Passed%9s Failed%9s Incomplete\n", "", "",
                         "", "", "");
    fla_test_output_info(
        "=====================================================================================\n");
    fla_test_output_info("%2sFLOAT%15s %8d%8s %8d%6s %8d%15s %d\n", "", "",
                         g_tests_passed[0] + g_tests_failed[0] + g_tests_incomplete[0], "",
                         g_tests_passed[0], "", g_tests_failed[0], "", g_tests_incomplete[0]);
    fla_test_output_info("%2sDOUBLE%14s %8d%8s %8d%6s %8d%15s %d\n", "", "",
                         g_tests_passed[1] + g_tests_failed[1] + g_tests_incomplete[0], "",
                         g_tests_passed[1], "", g_tests_failed[1], "", g_tests_incomplete[1]);
    fla_test_output_info("%2sCOMPLEX%13s %8d%8s %8d%6s %8d%15s %d\n", "", "",
                         g_tests_passed[2] + g_tests_failed[2] + g_tests_incomplete[0], "",
                         g_tests_passed[2], "", g_tests_failed[2], "", g_tests_incomplete[2]);
    fla_test_output_info("%2sDOUBLE COMPLEX%6s %8d%8s %8d%6s %8d%15s %d\n", "", "",
                         g_tests_passed[3] + g_tests_failed[3] + g_tests_incomplete[0], "",
                         g_tests_passed[3], "", g_tests_failed[3], "", g_tests_incomplete[3]);

    if(g_total_failed_tests > 0)
        printf("\n\nThere are failed tests, Please look at output log for more details\n\n");

    if(g_total_incomplete_tests > 0)
        printf("\n\nThere are some incomplete tests, Please look at the parameters passed to the "
               "API\n\n");
}

/*
 * This function dumps the runtimes to the file
 * specified in the dump_runtimes_file_name parameter
 * of the test_params_t structure.
 * The file is opened in append mode and the runtimes are
 * appended to the file.
 *
 * If the dump_runtimes_file_name is NULL, then
 * the function returns without doing anything.
 */
void fla_dump_runtimes_to_file(void *test_params, char *func_str, char datatype_char, integer p_cur, integer q_cur)
{
    test_params_t *params = (test_params_t *)test_params;
    test_runtime_ctx_t *ctx = &params->runtime_ctx;

    if(params->dump_runtimes_file_name == NULL)
    {
        return;
    }

    if(g_config_data == 1)
    {
        /* If config data is 1, then the test is invoked from the config file
           and the dump_runtimes_file_name is not used. */
        fla_test_output_error("Error: Dumping runtimes to file is not supported in config mode\n");
        return;
    }

    if(ctx->run_times_arr == NULL)
    {
        fla_test_output_error("Dump Error: Run times array is NULL\n");
        return;
    }

    FILE *fp = fopen(params->dump_runtimes_file_name, "a");
    if(fp == NULL)
    {
        fla_test_output_error("Error: Failed to open file %s for writing\n",
                              params->dump_runtimes_file_name);
        return;
    }

    fprintf(fp, "--------------------------------\n");
    fprintf(fp, "API: %s, Precision: %c, Size: %" FT_IS " x %" FT_IS ", Repeats: %zu\n", func_str, datatype_char, p_cur, q_cur, params->runtime_ctx.run_times_counter);
    fprintf(fp, "--------------------------------\n");
    for(integer i = 0; i < ctx->run_times_counter; ++i)
    {
        fprintf(fp, "%e\n", ctx->run_times_arr[i]);
    }
    fclose(fp);
    fla_test_output_info("Dumped runtimes to file %s\n", params->dump_runtimes_file_name);
}

/*
 * This function filters the runtimes based on the
 * filter_runtimes parameter of the test_params_t structure.
 * The function modifies the run_times_arr array in place.
 */
void fla_filter_runtimes(void *test_params)
{
    test_params_t *params = (test_params_t *)test_params;
    test_runtime_ctx_t *ctx = &params->runtime_ctx;

    if(params->outlier_multiplier == 0.0)
    {
        /* Don't do anything if filter is off */
        ctx->filtered_run_times_size = ctx->run_times_counter;
        return;
    }

    if(ctx->run_times_arr == NULL)
    {
        fla_test_output_error("Filter Error: Run times array is NULL. Filter not performed\n");
        ctx->filtered_run_times_size = ctx->run_times_counter;
        return;
    }

    double mean, stddev;

    mean = ctx->total_time / ctx->run_times_counter;
    get_stddev_of_array(DOUBLE, ctx->run_times_arr, &stddev, ctx->run_times_counter);

    double threshold = params->outlier_multiplier * stddev + mean;
    size_t insert_index = 0;

    for(size_t i = 0; i < ctx->run_times_counter; ++i)
    {
        if(ctx->run_times_arr[i] <= threshold)
        {
            ctx->run_times_arr[insert_index++] = ctx->run_times_arr[i];
        }
    }
    ctx->filtered_run_times_size = insert_index;
}