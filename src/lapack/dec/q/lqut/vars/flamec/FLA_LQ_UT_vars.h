/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_LQ_UT_unb_var1( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_LQ_UT_blk_var1( FLA_Obj A, FLA_Obj T, fla_lqut_t* cntl );
FLA_Error FLA_LQ_UT_opt_var1( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_LQ_UT_ops_var1( fla_dim_t m_A,
                              fla_dim_t n_A,
                              float* A, fla_dim_t rs_A, fla_dim_t cs_A,
                              float* t, fla_dim_t inc_t );
FLA_Error FLA_LQ_UT_opd_var1( fla_dim_t m_A,
                              fla_dim_t n_A,
                              double* A, fla_dim_t rs_A, fla_dim_t cs_A,
                              double* t, fla_dim_t inc_t );
FLA_Error FLA_LQ_UT_opc_var1( fla_dim_t m_A,
                              fla_dim_t n_A,
                              scomplex* A, fla_dim_t rs_A, fla_dim_t cs_A,
                              scomplex* t, fla_dim_t inc_t );
FLA_Error FLA_LQ_UT_opz_var1( fla_dim_t m_A,
                              fla_dim_t n_A,
                              dcomplex* A, fla_dim_t rs_A, fla_dim_t cs_A,
                              dcomplex* t, fla_dim_t inc_t );

FLA_Error FLA_LQ_UT_unb_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_LQ_UT_blk_var2( FLA_Obj A, FLA_Obj T, fla_lqut_t* cntl );
FLA_Error FLA_LQ_UT_opt_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_LQ_UT_ops_var2( fla_dim_t m_A,
                              fla_dim_t n_A,
                              float* A, fla_dim_t rs_A, fla_dim_t cs_A,
                              float* T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_LQ_UT_opd_var2( fla_dim_t m_A,
                              fla_dim_t n_A,
                              double* A, fla_dim_t rs_A, fla_dim_t cs_A,
                              double* T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_LQ_UT_opc_var2( fla_dim_t m_A,
                              fla_dim_t n_A,
                              scomplex* A, fla_dim_t rs_A, fla_dim_t cs_A,
                              scomplex* T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_LQ_UT_opz_var2( fla_dim_t m_A,
                              fla_dim_t n_A,
                              dcomplex* A, fla_dim_t rs_A, fla_dim_t cs_A,
                              dcomplex* T, fla_dim_t rs_T, fla_dim_t cs_T );

FLA_Error FLA_LQ_UT_blk_var3( FLA_Obj A, FLA_Obj T, fla_lqut_t* cntl );

