/******************************************************************************
 * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

#ifndef TEST_OUTPUT_ROUTINES_H
#define TEST_OUTPUT_ROUTINES_H

#include "test_common.h"

/* Print header in CLI mode if flag is true */
#define FLA_PRINT_HEADER_CLI                           \
    if(params->cli_print_header && g_config_data == 0) \
    {                                                  \
        fla_test_print_header(params);                 \
    }

void fla_test_output_info(char *message, ...);
void fla_test_output_error(char *message, ...);
void fla_test_print_header(void *test_params);
void fla_test_print_summary();
void fla_print_help(char *exe_name);
void fla_test_print_status(char *func_str, char datatype_char, integer sqr_inp, integer p_cur,
                           integer q_cur, double residual, double thresh, double perf,
                           void *params);
void fla_test_get_time_unit(char *scale, double measured_time, double *scale_factor,
                            integer time_unit);
integer fla_populate_stat_vals(double *stat_vals, void *test_params);

void fla_dump_runtimes_to_file(void *test_params, char *func_str, char datatype_char, integer p_cur, integer q_cur);

void fla_filter_runtimes(void *test_params);

#endif // TEST_OUTPUT_ROUTINES_H