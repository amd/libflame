/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_pivots_ln_opt_var1( FLA_Obj p, FLA_Obj A )
{
  FLA_Datatype datatype;
  fla_dim_t          n_A;
  fla_dim_t          rs_A, cs_A;
  fla_dim_t          inc_p;
  fla_dim_t          k1_0, k2_0;

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_p    = FLA_Obj_vector_inc( p );

  // Use zero-based indices.
  k1_0     = 0;
  k2_0     = ( fla_dim_t ) FLA_Obj_vector_dim( p ) - 1;

  switch ( datatype )
  {
    case FLA_INT:
    {
      fla_dim_t*   buff_A = FLA_INT_PTR( A );
      fla_dim_t*   buff_p = FLA_INT_PTR( p );

      FLA_Apply_pivots_ln_opi_var1( n_A,
                                    buff_A, rs_A, cs_A,
                                    k1_0,
                                    k2_0,
                                    buff_p, inc_p );

      break;
    }

    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      fla_dim_t*   buff_p = FLA_INT_PTR( p );

      FLA_Apply_pivots_ln_ops_var1( n_A,
                                    buff_A, rs_A, cs_A,
                                    k1_0,
                                    k2_0,
                                    buff_p, inc_p );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      fla_dim_t*    buff_p = FLA_INT_PTR( p );

      FLA_Apply_pivots_ln_opd_var1( n_A,
                                    buff_A, rs_A, cs_A,
                                    k1_0,
                                    k2_0,
                                    buff_p, inc_p );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      fla_dim_t*      buff_p = FLA_INT_PTR( p );

      FLA_Apply_pivots_ln_opc_var1( n_A,
                                    buff_A, rs_A, cs_A,
                                    k1_0,
                                    k2_0,
                                    buff_p, inc_p );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      fla_dim_t*      buff_p = FLA_INT_PTR( p );

      FLA_Apply_pivots_ln_opz_var1( n_A,
                                    buff_A, rs_A, cs_A,
                                    k1_0,
                                    k2_0,
                                    buff_p, inc_p );

      break;
    }
  }

  return FLA_SUCCESS;
}

FLA_Error FLA_Apply_pivots_ln_opi_var1( fla_dim_t n, 
                                        fla_dim_t* a, fla_dim_t a_rs, fla_dim_t a_cs, 
                                        fla_dim_t k1, 
                                        fla_dim_t k2, 
                                        fla_dim_t* p, fla_dim_t incp )
{
	fla_dim_t       temp;
	fla_dim_t*      a_i_0;
	fla_dim_t*      a_pi_0;
	fla_dim_t*      a_0_j;
	fla_dim_t*      a_i_j;
	fla_dim_t*      a_pi_j;
	fla_dim_t       i, j;
	fla_dim_t       i_begin, i_bound, i_inc;
	fla_dim_t       p_inc;

	// Handle both positive and negative increments for the pivot vector.
	if ( incp > 0 )
	{
		i_begin = k1;
		i_bound = k2 + 1;
		i_inc   = 1;
		p_inc   = 1*incp;
	}
	else // if ( incp < 0 )
	{
		i_begin = k2;
		i_bound = k1 - 1;
		i_inc   = -1;
		p_inc   = -1*incp;
	}

	// Optimize memory accesses depending on whether A is stored in
	// column-major or row-major order. That is, for column-major
	// matrices, we interchange all the elements in a single column
	// at a time. But for row-major matrices, we perform an entire
	// row interchange before moving to the next interchange. For
	// general storage, we decide based on which stride is closer
	// to one.
	if ( a_rs == 1 || a_rs < a_cs )
	{
		for ( j = 0; j < n; j++ )
		{
			a_0_j = a + j*a_cs;

			for ( i = i_begin; i != i_bound; i += i_inc )
			{
				a_i_j  = a_0_j + (              i )*a_rs;
				// Add i to shift from relative to absolute index.
				a_pi_j = a_0_j + ( p[i*p_inc] + i )*a_rs;

				temp    = *a_pi_j;
				*a_pi_j = *a_i_j;
				*a_i_j  = temp;
			}
		}
	}
	else // if ( a_cs == 1 || a_cs < a_rs )
	{
		for ( i = i_begin; i != i_bound; i += i_inc )
		{
			a_i_0  = a + (              i )*a_rs;
			// Add i to shift from relative to absolute index.
			a_pi_0 = a + ( p[i*p_inc] + i )*a_rs;

			for ( j = 0; j < n; j++ )
			{
				a_i_j  = a_i_0 + j*a_cs;
				a_pi_j = a_pi_0 + j*a_cs;

				temp    = *a_pi_j;
				*a_pi_j = *a_i_j;
				*a_i_j  = temp;
			}
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Apply_pivots_ln_ops_var1( fla_dim_t n, 
                                        float* a, fla_dim_t a_rs, fla_dim_t a_cs, 
                                        fla_dim_t k1, 
                                        fla_dim_t k2, 
                                        fla_dim_t* p, fla_dim_t incp )
{
	float     temp;
	float*    a_i_0;
	float*    a_pi_0;
	float*    a_0_j;
	float*    a_i_j;
	float*    a_pi_j;
	fla_dim_t       i, j;
	fla_dim_t       i_begin, i_bound, i_inc;
	fla_dim_t       p_inc;

	// Handle both positive and negative increments for the pivot vector.
	if ( incp > 0 )
	{
		i_begin = k1;
		i_bound = k2 + 1;
		i_inc   = 1;
		p_inc   = 1*incp;
	}
	else // if ( incp < 0 )
	{
		i_begin = k2;
		i_bound = k1 - 1;
		i_inc   = -1;
		p_inc   = -1*incp;
	}

	// Optimize memory accesses depending on whether A is stored in
	// column-major or row-major order. That is, for column-major
	// matrices, we interchange all the elements in a single column
	// at a time. But for row-major matrices, we perform an entire
	// row interchange before moving to the next interchange. For
	// general storage, we decide based on which stride is closer
	// to one.
	if ( a_rs == 1 || a_rs < a_cs )
	{
		for ( j = 0; j < n; j++ )
		{
			a_0_j = a + j*a_cs;

			for ( i = i_begin; i != i_bound; i += i_inc )
			{
				a_i_j  = a_0_j + (              i )*a_rs;
				// Add i to shift from relative to absolute index.
				a_pi_j = a_0_j + ( p[i*p_inc] + i )*a_rs;

				temp    = *a_pi_j;
				*a_pi_j = *a_i_j;
				*a_i_j  = temp;
			}
		}
	}
	else // if ( a_cs == 1 || a_cs < a_rs )
	{
		for ( i = i_begin; i != i_bound; i += i_inc )
		{
			a_i_0  = a + (              i )*a_rs;
			// Add i to shift from relative to absolute index.
			a_pi_0 = a + ( p[i*p_inc] + i )*a_rs;

			for ( j = 0; j < n; j++ )
			{
				a_i_j  = a_i_0 + j*a_cs;
				a_pi_j = a_pi_0 + j*a_cs;

				temp    = *a_pi_j;
				*a_pi_j = *a_i_j;
				*a_i_j  = temp;
			}
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Apply_pivots_ln_opd_var1( fla_dim_t n,
                                        double* a, fla_dim_t a_rs, fla_dim_t a_cs,
                                        fla_dim_t k1,
                                        fla_dim_t k2,
                                        fla_dim_t* p, fla_dim_t incp )
{
	double    temp;
	double*   a_i_0;
	double*   a_pi_0;
	double*   a_0_j;
	double*   a_i_j;
	double*   a_pi_j;
	fla_dim_t       i, j;
	fla_dim_t       i_begin, i_bound, i_inc;
	fla_dim_t       p_inc;

	// Handle both positive and negative increments for the pivot vector.
	if ( incp > 0 )
	{
		i_begin = k1;
		i_bound = k2 + 1;
		i_inc   = 1;
		p_inc   = 1*incp;
	}
	else // if ( incp < 0 )
	{
		i_begin = k2;
		i_bound = k1 - 1;
		i_inc   = -1;
		p_inc   = -1*incp;
	}

	// Optimize memory accesses depending on whether A is stored in
	// column-major or row-major order. That is, for column-major
	// matrices, we interchange all the elements in a single column
	// at a time. But for row-major matrices, we perform an entire
	// row interchange before moving to the next interchange. For
	// general storage, we decide based on which stride is closer
	// to one.
	if ( a_rs == 1 || a_rs < a_cs )
	{
		for ( j = 0; j < n; j++ )
		{
			a_0_j = a + j*a_cs;

			for ( i = i_begin; i != i_bound; i += i_inc )
			{
				a_i_j  = a_0_j + (              i )*a_rs;
				// Add i to shift from relative to absolute index.
				a_pi_j = a_0_j + ( p[i*p_inc] + i )*a_rs;

				temp    = *a_pi_j;
				*a_pi_j = *a_i_j;
				*a_i_j  = temp;
			}
		}
	}
	else // if ( a_cs == 1 || a_cs < a_rs )
	{
		for ( i = i_begin; i != i_bound; i += i_inc )
		{
			a_i_0  = a + (              i )*a_rs;
			// Add i to shift from relative to absolute index.
			a_pi_0 = a + ( p[i*p_inc] + i )*a_rs;

			for ( j = 0; j < n; j++ )
			{
				a_i_j  = a_i_0 + j*a_cs;
				a_pi_j = a_pi_0 + j*a_cs;

				temp    = *a_pi_j;
				*a_pi_j = *a_i_j;
				*a_i_j  = temp;
			}
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Apply_pivots_ln_opc_var1( fla_dim_t n,
                                        scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs,
                                        fla_dim_t k1,
                                        fla_dim_t k2,
                                        fla_dim_t* p, fla_dim_t incp )
{
	scomplex  temp;
	scomplex* a_i_0;
	scomplex* a_pi_0;
	scomplex* a_0_j;
	scomplex* a_i_j;
	scomplex* a_pi_j;
	fla_dim_t       i, j;
	fla_dim_t       i_begin, i_bound, i_inc;
	fla_dim_t       p_inc;

	// Handle both positive and negative increments for the pivot vector.
	if ( incp > 0 )
	{
		i_begin = k1;
		i_bound = k2 + 1;
		i_inc   = 1;
		p_inc   = 1*incp;
	}
	else // if ( incp < 0 )
	{
		i_begin = k2;
		i_bound = k1 - 1;
		i_inc   = -1;
		p_inc   = -1*incp;
	}

	// Optimize memory accesses depending on whether A is stored in
	// column-major or row-major order. That is, for column-major
	// matrices, we interchange all the elements in a single column
	// at a time. But for row-major matrices, we perform an entire
	// row interchange before moving to the next interchange. For
	// general storage, we decide based on which stride is closer
	// to one.
	if ( a_rs == 1 || a_rs < a_cs )
	{
		for ( j = 0; j < n; j++ )
		{
			a_0_j = a + j*a_cs;

			for ( i = i_begin; i != i_bound; i += i_inc )
			{
				a_i_j  = a_0_j + (              i )*a_rs;
				// Add i to shift from relative to absolute index.
				a_pi_j = a_0_j + ( p[i*p_inc] + i )*a_rs;

				temp    = *a_pi_j;
				*a_pi_j = *a_i_j;
				*a_i_j  = temp;
			}
		}
	}
	else // if ( a_cs == 1 || a_cs < a_rs )
	{
		for ( i = i_begin; i != i_bound; i += i_inc )
		{
			a_i_0  = a + (              i )*a_rs;
			// Add i to shift from relative to absolute index.
			a_pi_0 = a + ( p[i*p_inc] + i )*a_rs;

			for ( j = 0; j < n; j++ )
			{
				a_i_j  = a_i_0 + j*a_cs;
				a_pi_j = a_pi_0 + j*a_cs;

				temp    = *a_pi_j;
				*a_pi_j = *a_i_j;
				*a_i_j  = temp;
			}
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Apply_pivots_ln_opz_var1( fla_dim_t n,
                                        dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs,
                                        fla_dim_t k1,
                                        fla_dim_t k2,
                                        fla_dim_t* p, fla_dim_t incp )
{
	dcomplex  temp;
	dcomplex* a_i_0;
	dcomplex* a_pi_0;
	dcomplex* a_0_j;
	dcomplex* a_i_j;
	dcomplex* a_pi_j;
	fla_dim_t       i, j;
	fla_dim_t       i_begin, i_bound, i_inc;
	fla_dim_t       p_inc;

	// Handle both positive and negative increments for the pivot vector.
	if ( incp > 0 )
	{
		i_begin = k1;
		i_bound = k2 + 1;
		i_inc   = 1;
		p_inc   = 1*incp;
	}
	else // if ( incp < 0 )
	{
		i_begin = k2;
		i_bound = k1 - 1;
		i_inc   = -1;
		p_inc   = -1*incp;
	}

	// Optimize memory accesses depending on whether A is stored in
	// column-major or row-major order. That is, for column-major
	// matrices, we interchange all the elements in a single column
	// at a time. But for row-major matrices, we perform an entire
	// row interchange before moving to the next interchange. For
	// general storage, we decide based on which stride is closer
	// to one.
	if ( a_rs == 1 || a_rs < a_cs )
	{
		for ( j = 0; j < n; j++ )
		{
			a_0_j = a + j*a_cs;

			for ( i = i_begin; i != i_bound; i += i_inc )
			{
				a_i_j  = a_0_j + (              i )*a_rs;
				// Add i to shift from relative to absolute index.
				a_pi_j = a_0_j + ( p[i*p_inc] + i )*a_rs;

				temp    = *a_pi_j;
				*a_pi_j = *a_i_j;
				*a_i_j  = temp;
			}
		}
	}
	else // if ( a_cs == 1 || a_cs < a_rs )
	{
		for ( i = i_begin; i != i_bound; i += i_inc )
		{
			a_i_0  = a + (              i )*a_rs;
			// Add i to shift from relative to absolute index.
			a_pi_0 = a + ( p[i*p_inc] + i )*a_rs;

			for ( j = 0; j < n; j++ )
			{
				a_i_j  = a_i_0 + j*a_cs;
				a_pi_j = a_pi_0 + j*a_cs;

				temp    = *a_pi_j;
				*a_pi_j = *a_i_j;
				*a_i_j  = temp;
			}
		}
	}

	return FLA_SUCCESS;
}

