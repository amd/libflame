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

int bl1_vector_dim( fla_dim_t m, fla_dim_t n )
{
	if ( m == 1 ) return n;
	else          return m;
}

int bl1_vector_inc( trans1_t trans, fla_dim_t m, fla_dim_t n, fla_dim_t rs, fla_dim_t cs )
{
	if ( bl1_does_notrans( trans ) )
	{
		if ( m == 1 ) return cs;
		else          return rs;
	}
	else // if ( bl1_does_trans( trans ) )
	{
		if ( m == 1 ) return rs;
		else          return cs;
	}
}
