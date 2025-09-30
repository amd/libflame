/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
*     Modifications Copyright (c) 2023 Advanced Micro Devices, Inc.  All rights reserved.
*/
#include "blis1.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif

// --- storage-related ---------------------------------------------------------

void bl1_check_storage_3m( fla_dim_t a_rs, fla_dim_t a_cs, fla_dim_t b_rs, fla_dim_t b_cs, fla_dim_t c_rs, fla_dim_t c_cs )
{
	if ( bl1_is_gen_storage( a_rs, a_cs ) ||
	     bl1_is_gen_storage( b_rs, b_cs ) ||
	     bl1_is_gen_storage( c_rs, c_cs ) )
	{
		bl1_abort_msg( "Function or conditional branch/case not yet implemented." );
	}
}

void bl1_check_storage_2m( fla_dim_t a_rs, fla_dim_t a_cs, fla_dim_t b_rs, fla_dim_t b_cs )
{
	if ( bl1_is_gen_storage( a_rs, a_cs ) ||
	     bl1_is_gen_storage( b_rs, b_cs ) )
	{
		bl1_abort_msg( "Function or conditional branch/case not yet implemented." );
	}
}

