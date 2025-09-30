/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLA_Tridiag_UT_l_blk_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_unb_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_unb_var1( FLA_Obj A, FLA_Obj T );

FLA_Error FLA_Tridiag_UT_l_blk_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_blf_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_unb_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_unb_var2( FLA_Obj A, FLA_Obj T );

FLA_Error FLA_Tridiag_UT_l_blk_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_blf_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_unb_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_unb_var3( FLA_Obj A, FLA_Obj Z, FLA_Obj T );

FLA_Error FLA_Tridiag_UT_l_opt_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_opt_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ops_var1( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          float* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opd_var1( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          double* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opc_var1( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          scomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opz_var1( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          dcomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );

FLA_Error FLA_Tridiag_UT_l_opt_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_opt_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ops_var2( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          float* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opd_var2( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          double* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opc_var2( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          scomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opz_var2( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          dcomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );

FLA_Error FLA_Tridiag_UT_l_opt_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_opt_var3( FLA_Obj A, FLA_Obj Z, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ops_var3( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          float* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                          float* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opd_var3( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          double* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                          double* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opc_var3( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          scomplex* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                          scomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opz_var3( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          dcomplex* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                          dcomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );

FLA_Error FLA_Tridiag_UT_l_ofu_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ofu_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ofs_var1( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          float* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofd_var1( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          double* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofc_var1( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          scomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofz_var1( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          dcomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );

FLA_Error FLA_Tridiag_UT_l_ofu_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ofu_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ofs_var2( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          float* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofd_var2( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          double* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofc_var2( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          scomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofz_var2( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          dcomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );

FLA_Error FLA_Tridiag_UT_l_ofu_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ofu_var3( FLA_Obj A, FLA_Obj Z, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ofs_var3( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          float* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z,
                                          float* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofd_var3( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          double* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z,
                                          double* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofc_var3( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          scomplex* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z,
                                          scomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofz_var3( fla_dim_t m_A,
                                          fla_dim_t m_T,
                                          dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          dcomplex* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z,
                                          dcomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T );

// --- Fused operations ---

FLA_Error FLA_Fused_Her2_Ax_l_opt_var1( FLA_Obj alpha, FLA_Obj u, FLA_Obj z, FLA_Obj A, FLA_Obj x, FLA_Obj w );
FLA_Error FLA_Fused_Her2_Ax_l_ops_var1( fla_dim_t m_A,
                                        float* buff_alpha, 
                                        float* buff_u, fla_dim_t inc_u, 
                                        float* buff_z, fla_dim_t inc_z, 
                                        float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                        float* buff_x, fla_dim_t inc_x, 
                                        float* buff_w, fla_dim_t inc_w );
FLA_Error FLA_Fused_Her2_Ax_l_opd_var1( fla_dim_t m_A,
                                        double* buff_alpha, 
                                        double* buff_u, fla_dim_t inc_u, 
                                        double* buff_z, fla_dim_t inc_z, 
                                        double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                        double* buff_x, fla_dim_t inc_x, 
                                        double* buff_w, fla_dim_t inc_w );
FLA_Error FLA_Fused_Her2_Ax_l_opc_var1( fla_dim_t m_A,
                                        scomplex* buff_alpha, 
                                        scomplex* buff_u, fla_dim_t inc_u, 
                                        scomplex* buff_z, fla_dim_t inc_z, 
                                        scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                        scomplex* buff_x, fla_dim_t inc_x, 
                                        scomplex* buff_w, fla_dim_t inc_w );
FLA_Error FLA_Fused_Her2_Ax_l_opz_var1( fla_dim_t m_A,
                                        dcomplex* buff_alpha, 
                                        dcomplex* buff_u, fla_dim_t inc_u, 
                                        dcomplex* buff_z, fla_dim_t inc_z, 
                                        dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                        dcomplex* buff_x, fla_dim_t inc_x, 
                                        dcomplex* buff_w, fla_dim_t inc_w );

FLA_Error FLA_Fused_UZhu_ZUhu_opt_var1( FLA_Obj delta, FLA_Obj U, FLA_Obj Z, FLA_Obj t, FLA_Obj u, FLA_Obj w );
FLA_Error FLA_Fused_UZhu_ZUhu_ops_var1( fla_dim_t m_U,
                                        fla_dim_t n_U,
                                        float* buff_delta, 
                                        float* buff_U, fla_dim_t rs_U, fla_dim_t cs_U, 
                                        float* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                        float* buff_t, fla_dim_t inc_t, 
                                        float* buff_u, fla_dim_t inc_u, 
                                        float* buff_w, fla_dim_t inc_w );
FLA_Error FLA_Fused_UZhu_ZUhu_opd_var1( fla_dim_t m_U,
                                        fla_dim_t n_U,
                                        double* buff_delta, 
                                        double* buff_U, fla_dim_t rs_U, fla_dim_t cs_U, 
                                        double* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                        double* buff_t, fla_dim_t inc_t, 
                                        double* buff_u, fla_dim_t inc_u, 
                                        double* buff_w, fla_dim_t inc_w );
FLA_Error FLA_Fused_UZhu_ZUhu_opc_var1( fla_dim_t m_U,
                                        fla_dim_t n_U,
                                        scomplex* buff_delta, 
                                        scomplex* buff_U, fla_dim_t rs_U, fla_dim_t cs_U, 
                                        scomplex* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                        scomplex* buff_t, fla_dim_t inc_t, 
                                        scomplex* buff_u, fla_dim_t inc_u, 
                                        scomplex* buff_w, fla_dim_t inc_w );
FLA_Error FLA_Fused_UZhu_ZUhu_opz_var1( fla_dim_t m_U,
                                        fla_dim_t n_U,
                                        dcomplex* buff_delta, 
                                        dcomplex* buff_U, fla_dim_t rs_U, fla_dim_t cs_U, 
                                        dcomplex* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                        dcomplex* buff_t, fla_dim_t inc_t, 
                                        dcomplex* buff_u, fla_dim_t inc_u, 
                                        dcomplex* buff_w, fla_dim_t inc_w );
