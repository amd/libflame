/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#ifdef FLA_OPENMP_MULTITHREADING

#ifndef API_ID_DEFINED
#define API_ID_DEFINED

#include "FLA_type_defs.h"

/* API ID */
typedef enum
{
    FLA_LABRD = 0,
    FLA_ORMQR,
    FLA_ORMLQ
} API_ID;
#endif


void FLA_Thread_get_subrange( int thread_ID, int num_threads, fla_dim_t range, fla_dim_t *sub_range, fla_dim_t *index );
void FLA_Thread_optimum( API_ID family, int *actual_num_threads);

#endif