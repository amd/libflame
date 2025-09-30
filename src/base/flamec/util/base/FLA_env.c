/*
    Copyright (c) 2021-2023 Advanced Micro Devices, Inc. All rights reserved.
    Mar 16, 2021
*/

#include "FLAME.h"

fla_dim_t FLA_env_get_var( const char* env, fla_dim_t fallback )
{
	fla_dim_t r_val;
	char*  str;

	// Query the environment variable and store the result in str.
	str = getenv( env );

	// Set the return value based on the string obtained from getenv().
	if ( str != NULL )
	{
		// If there was no error, convert the string to an fla_dim_t and
		// prepare to return that fla_dim_t.
		r_val = ( fla_dim_t )strtol( str, NULL, 10 );
	}
	else
	{
		// If there was an error, use the "fallback" as the return value.
		r_val = fallback;
	}

	return r_val;
}

fla_dim_t FLASH_get_num_threads( fla_dim_t fallback )
{

    fla_dim_t fla_threads;
    extern int fla_thread_get_num_threads(void);

    fla_threads = fla_thread_get_num_threads();

    if( fla_threads <= 0 )
    {
        fla_threads = fallback;
    }

    return fla_threads;
}

