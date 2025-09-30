/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLA_UDdate_UT_blk_var1( FLA_Obj R,
                                  FLA_Obj C,
                                  FLA_Obj D, FLA_Obj T, fla_uddateut_t* cntl );
FLA_Error FLA_UDdate_UT_blk_var2( FLA_Obj R,
                                  FLA_Obj C,
                                  FLA_Obj D, FLA_Obj T, fla_uddateut_t* cntl );

FLA_Error FLA_UDdate_UT_unb_var1( FLA_Obj R,
                                  FLA_Obj C,
                                  FLA_Obj D, FLA_Obj T );

FLA_Error FLA_UDdate_UT_opt_var1( FLA_Obj R,
                                  FLA_Obj C,
                                  FLA_Obj D, FLA_Obj T );
FLA_Error FLA_UDdate_UT_ops_var1( fla_dim_t mn_RT,
                                  fla_dim_t m_C,
                                  fla_dim_t m_D,
                                  float* R, fla_dim_t rs_R, fla_dim_t cs_R,
                                  float* C, fla_dim_t rs_C, fla_dim_t cs_C,
                                  float* D, fla_dim_t rs_D, fla_dim_t cs_D,
                                  float* T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_UDdate_UT_opd_var1( fla_dim_t mn_RT,
                                  fla_dim_t m_C,
                                  fla_dim_t m_D,
                                  double* R, fla_dim_t rs_R, fla_dim_t cs_R,
                                  double* C, fla_dim_t rs_C, fla_dim_t cs_C,
                                  double* D, fla_dim_t rs_D, fla_dim_t cs_D,
                                  double* T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_UDdate_UT_opc_var1( fla_dim_t mn_RT,
                                  fla_dim_t m_C,
                                  fla_dim_t m_D,
                                  scomplex* R, fla_dim_t rs_R, fla_dim_t cs_R,
                                  scomplex* C, fla_dim_t rs_C, fla_dim_t cs_C,
                                  scomplex* D, fla_dim_t rs_D, fla_dim_t cs_D,
                                  scomplex* T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_UDdate_UT_opz_var1( fla_dim_t mn_RT,
                                  fla_dim_t m_C,
                                  fla_dim_t m_D,
                                  dcomplex* R, fla_dim_t rs_R, fla_dim_t cs_R,
                                  dcomplex* C, fla_dim_t rs_C, fla_dim_t cs_C,
                                  dcomplex* D, fla_dim_t rs_D, fla_dim_t cs_D,
                                  dcomplex* T, fla_dim_t rs_T, fla_dim_t cs_T );

