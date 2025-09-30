/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Accum_T_UT_fc_unb_var1( FLA_Obj A, FLA_Obj t, FLA_Obj T );
FLA_Error FLA_Accum_T_UT_fc_blk_var2( FLA_Obj A, FLA_Obj t, FLA_Obj T );

FLA_Error FLA_Accum_T_UT_fc_opt_var1( FLA_Obj A, FLA_Obj t, FLA_Obj T );

FLA_Error FLA_Accum_T_UT_fc_ops_var1( fla_dim_t m_A,
                                      fla_dim_t n_AT,
                                      float* A, fla_dim_t rs_A, fla_dim_t cs_A,
                                      fla_dim_t m_t, 
                                      float* t, fla_dim_t inc_t,
                                      float* T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Accum_T_UT_fc_opd_var1( fla_dim_t m_A,
                                      fla_dim_t n_AT,
                                      double* A, fla_dim_t rs_A, fla_dim_t cs_A,
                                      fla_dim_t m_t, 
                                      double* t, fla_dim_t inc_t,
                                      double* T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Accum_T_UT_fc_opc_var1( fla_dim_t m_A,
                                      fla_dim_t n_AT,
                                      scomplex* A, fla_dim_t rs_A, fla_dim_t cs_A,
                                      fla_dim_t m_t, 
                                      scomplex* t, fla_dim_t inc_t,
                                      scomplex* T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Accum_T_UT_fc_opz_var1( fla_dim_t m_A,
                                      fla_dim_t n_AT,
                                      dcomplex* A, fla_dim_t rs_A, fla_dim_t cs_A,
                                      fla_dim_t m_t, 
                                      dcomplex* t, fla_dim_t inc_t,
                                      dcomplex* T, fla_dim_t rs_T, fla_dim_t cs_T );
