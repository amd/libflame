/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_G_rf_blk_var9( FLA_Obj G, FLA_Obj A, fla_dim_t b_alg )
{
	FLA_Datatype datatype;
	fla_dim_t          k_G, m_A, n_A;
	fla_dim_t          rs_G, cs_G;
	fla_dim_t          rs_A, cs_A;

	datatype = FLA_Obj_datatype( A );

	k_G      = FLA_Obj_width( G );
	m_A      = FLA_Obj_length( A );
	n_A      = FLA_Obj_width( A );

	rs_G     = FLA_Obj_row_stride( G );
	cs_G     = FLA_Obj_col_stride( G );

	rs_A     = FLA_Obj_row_stride( A );
	cs_A     = FLA_Obj_col_stride( A );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			scomplex* buff_G = ( scomplex* ) FLA_COMPLEX_PTR( G );
			float*    buff_A = ( float*    ) FLA_FLOAT_PTR( A );

			FLA_Apply_G_rf_bls_var9( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A,
			                         b_alg );

			break;
		}

		case FLA_DOUBLE:
		{
			dcomplex* buff_G = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( G );
			double*   buff_A = ( double*   ) FLA_DOUBLE_PTR( A );

			FLA_Apply_G_rf_bld_var9( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A,
			                         b_alg );

			break;
		}

		case FLA_COMPLEX:
		{
			scomplex* buff_G = ( scomplex* ) FLA_COMPLEX_PTR( G );
			scomplex* buff_A = ( scomplex* ) FLA_COMPLEX_PTR( A );

			FLA_Apply_G_rf_blc_var9( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A,
			                         b_alg );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			dcomplex* buff_G = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( G );
			dcomplex* buff_A = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );

			FLA_Apply_G_rf_blz_var9( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A,
			                         b_alg );

			break;
		}
	}

	return FLA_SUCCESS;
}


FLA_Error FLA_Apply_G_rf_bls_var9( fla_dim_t       k_G,
                                   fla_dim_t       m_A,
                                   fla_dim_t       n_A,
                                   scomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                                   float*    buff_A, fla_dim_t rs_A, fla_dim_t cs_A,
                                   fla_dim_t       b_alg )
{
	fla_dim_t i;
	fla_dim_t b = 0;

	for ( i = 0; i < m_A; i += b )
	{
		float*    A1      = buff_A + (0  )*cs_A + (i  )*rs_A;
		fla_dim_t       m_ahead = fla_max( 0, m_A - i );

		b = fla_min( b_alg, m_ahead );

		//FLA_Apply_G_rf_ops_var9( k_G,
		FLA_Apply_G_rf_ass_var9( k_G,
		                         b,
		                         n_A,
		                         buff_G, rs_G, cs_G,
		                         A1,     rs_A, cs_A );
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_bld_var9( fla_dim_t       k_G,
                                   fla_dim_t       m_A,
                                   fla_dim_t       n_A,
                                   dcomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                                   double*   buff_A, fla_dim_t rs_A, fla_dim_t cs_A,
                                   fla_dim_t       b_alg )
{
	fla_dim_t i;
	fla_dim_t b = 0;

	for ( i = 0; i < m_A; i += b )
	{
		double*   A1      = buff_A + (0  )*cs_A + (i  )*rs_A;
		fla_dim_t       m_ahead = fla_max( 0, m_A - i );

		b = fla_min( b_alg, m_ahead );

		//FLA_Apply_G_rf_opd_var9( k_G,
		FLA_Apply_G_rf_asd_var9( k_G,
		                         b,
		                         n_A,
		                         buff_G, rs_G, cs_G,
		                         A1,     rs_A, cs_A );
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_blc_var9( fla_dim_t       k_G,
                                   fla_dim_t       m_A,
                                   fla_dim_t       n_A,
                                   scomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                                   scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A,
                                   fla_dim_t       b_alg )
{
	fla_dim_t i;
	fla_dim_t b = 0;

	for ( i = 0; i < m_A; i += b )
	{
		scomplex* A1      = buff_A + (0  )*cs_A + (i  )*rs_A;
		fla_dim_t       m_ahead = fla_max( 0, m_A - i );

		b = fla_min( b_alg, m_ahead );

		//FLA_Apply_G_rf_opc_var9( k_G,
		FLA_Apply_G_rf_asc_var9( k_G,
		                         b,
		                         n_A,
		                         buff_G, rs_G, cs_G,
		                         A1,     rs_A, cs_A );
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_blz_var9( fla_dim_t       k_G,
                                   fla_dim_t       m_A,
                                   fla_dim_t       n_A,
                                   dcomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                                   dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A,
                                   fla_dim_t       b_alg )
{
	fla_dim_t i;
	fla_dim_t b = 0;

	for ( i = 0; i < m_A; i += b )
	{
		dcomplex* A1      = buff_A + (0  )*cs_A + (i  )*rs_A;
		fla_dim_t       m_ahead = fla_max( 0, m_A - i );

		b = fla_min( b_alg, m_ahead );

		//FLA_Apply_G_rf_opz_var9( k_G,
		FLA_Apply_G_rf_asz_var9( k_G,
		                         b,
		                         n_A,
		                         buff_G, rs_G, cs_G,
		                         A1,     rs_A, cs_A );
	}

	return FLA_SUCCESS;
}
