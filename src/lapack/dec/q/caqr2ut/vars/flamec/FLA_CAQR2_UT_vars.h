/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLA_CAQR2_UT_blk_var1( FLA_Obj U,
                                 FLA_Obj D, FLA_Obj T, fla_caqr2ut_t* cntl );
FLA_Error FLA_CAQR2_UT_blk_var2( FLA_Obj U,
                                 FLA_Obj D, FLA_Obj T, fla_caqr2ut_t* cntl );

FLA_Error FLA_CAQR2_UT_unb_var1( FLA_Obj U,
                                 FLA_Obj D, FLA_Obj T );

FLA_Error FLA_CAQR2_UT_opt_var1( FLA_Obj U,
                                 FLA_Obj D, FLA_Obj T );
FLA_Error FLA_CAQR2_UT_ops_var1( fla_dim_t m_UT,
                                 fla_dim_t m_D,
                                 float* U, fla_dim_t rs_U, fla_dim_t cs_U,
                                 float* D, fla_dim_t rs_D, fla_dim_t cs_D,
                                 float* T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_CAQR2_UT_opd_var1( fla_dim_t m_UT,
                                 fla_dim_t m_D,
                                 double* U, fla_dim_t rs_U, fla_dim_t cs_U,
                                 double* D, fla_dim_t rs_D, fla_dim_t cs_D,
                                 double* T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_CAQR2_UT_opc_var1( fla_dim_t m_UT,
                                 fla_dim_t m_D,
                                 scomplex* U, fla_dim_t rs_U, fla_dim_t cs_U,
                                 scomplex* D, fla_dim_t rs_D, fla_dim_t cs_D,
                                 scomplex* T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_CAQR2_UT_opz_var1( fla_dim_t m_UT,
                                 fla_dim_t m_D,
                                 dcomplex* U, fla_dim_t rs_U, fla_dim_t cs_U,
                                 dcomplex* D, fla_dim_t rs_D, fla_dim_t cs_D,
                                 dcomplex* T, fla_dim_t rs_T, fla_dim_t cs_T );

