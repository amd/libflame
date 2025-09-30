/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLA_Bidiag_UT_u_unb_var1( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blk_var1( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_step_unb_var1( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );

FLA_Error FLA_Bidiag_UT_u_unb_var2( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blk_var2( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blf_var2( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_step_unb_var2( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );

FLA_Error FLA_Bidiag_UT_u_unb_var3( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blk_var3( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blf_var3( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_step_unb_var3( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );

FLA_Error FLA_Bidiag_UT_u_unb_var4( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blk_var4( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blf_var4( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_step_unb_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj TU, FLA_Obj TV );

FLA_Error FLA_Bidiag_UT_u_unb_var5( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blk_var5( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_step_unb_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj TU, FLA_Obj TV );

FLA_Error FLA_Bidiag_UT_u_opt_var1( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_opt_var1( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ops_var1( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         float* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         float* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opd_var1( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         double* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         double* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opc_var1( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         scomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         scomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opz_var1( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         dcomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         dcomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );

FLA_Error FLA_Bidiag_UT_u_opt_var2( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_opt_var2( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ops_var2( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         float* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         float* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opd_var2( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         double* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         double* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opc_var2( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         scomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         scomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opz_var2( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         dcomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         dcomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );

FLA_Error FLA_Bidiag_UT_u_opt_var3( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_opt_var3( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ops_var3( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         float* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         float* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opd_var3( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         double* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         double* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opc_var3( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         scomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         scomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opz_var3( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         dcomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         dcomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );

FLA_Error FLA_Bidiag_UT_u_opt_var4( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_opt_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ops_var4( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         float* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                         float* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                         float* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         float* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opd_var4( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         double* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                         double* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                         double* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         double* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opc_var4( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         scomplex* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                         scomplex* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                         scomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         scomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opz_var4( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         dcomplex* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                         dcomplex* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                         dcomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         dcomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );

FLA_Error FLA_Bidiag_UT_u_opt_var5( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_opt_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ops_var5( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         float* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                         float* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                         float* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         float* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opd_var5( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         double* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                         double* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                         double* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         double* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opc_var5( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         scomplex* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                         scomplex* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                         scomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         scomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opz_var5( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         dcomplex* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                         dcomplex* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                         dcomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         dcomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );


FLA_Error FLA_Bidiag_UT_u_ofu_var2( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ofu_var2( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ofs_var2( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         float* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         float* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofd_var2( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         double* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         double* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofc_var2( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         scomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         scomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofz_var2( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         dcomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         dcomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );

FLA_Error FLA_Bidiag_UT_u_ofu_var3( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ofu_var3( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ofs_var3( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         float* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         float* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofd_var3( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         double* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         double* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofc_var3( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         scomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         scomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofz_var3( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         dcomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         dcomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );

FLA_Error FLA_Bidiag_UT_u_ofu_var4( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ofu_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ofs_var4( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         float* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                         float* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                         float* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         float* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofd_var4( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         double* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                         double* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                         double* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         double* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofc_var4( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         scomplex* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                         scomplex* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                         scomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         scomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofz_var4( fla_dim_t m_A,
                                         fla_dim_t n_A,
                                         fla_dim_t m_TS,
                                         dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                         dcomplex* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                         dcomplex* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                         dcomplex* buff_T, fla_dim_t rs_T, fla_dim_t cs_T, 
                                         dcomplex* buff_S, fla_dim_t rs_S, fla_dim_t cs_S );

// --- Fused operations ---

FLA_Error FLA_Fused_Gerc2_opt_var1( FLA_Obj alpha, FLA_Obj u, FLA_Obj y, FLA_Obj z, FLA_Obj v, FLA_Obj A );
FLA_Error FLA_Fused_Gerc2_ops_var1( fla_dim_t m_A,
                                    fla_dim_t n_A,
                                    float* buff_alpha, 
                                    float* buff_u, fla_dim_t inc_u, 
                                    float* buff_y, fla_dim_t inc_y, 
                                    float* buff_z, fla_dim_t inc_z, 
                                    float* buff_v, fla_dim_t inc_v, 
                                    float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A ); 
FLA_Error FLA_Fused_Gerc2_opd_var1( fla_dim_t m_A,
                                    fla_dim_t n_A,
                                    double* buff_alpha, 
                                    double* buff_u, fla_dim_t inc_u, 
                                    double* buff_y, fla_dim_t inc_y, 
                                    double* buff_z, fla_dim_t inc_z, 
                                    double* buff_v, fla_dim_t inc_v, 
                                    double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A ); 
FLA_Error FLA_Fused_Gerc2_opc_var1( fla_dim_t m_A,
                                    fla_dim_t n_A,
                                    scomplex* buff_alpha, 
                                    scomplex* buff_u, fla_dim_t inc_u, 
                                    scomplex* buff_y, fla_dim_t inc_y, 
                                    scomplex* buff_z, fla_dim_t inc_z, 
                                    scomplex* buff_v, fla_dim_t inc_v, 
                                    scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A ); 
FLA_Error FLA_Fused_Gerc2_opz_var1( fla_dim_t m_A,
                                    fla_dim_t n_A,
                                    dcomplex* buff_alpha, 
                                    dcomplex* buff_u, fla_dim_t inc_u, 
                                    dcomplex* buff_y, fla_dim_t inc_y, 
                                    dcomplex* buff_z, fla_dim_t inc_z, 
                                    dcomplex* buff_v, fla_dim_t inc_v, 
                                    dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A ); 


FLA_Error FLA_Fused_Ahx_Axpy_Ax_opt_var1( FLA_Obj A, FLA_Obj u, FLA_Obj tau, FLA_Obj a, FLA_Obj beta, FLA_Obj y, FLA_Obj w );
FLA_Error FLA_Fused_Ahx_Axpy_Ax_ops_var1( fla_dim_t m_A,
                                          fla_dim_t n_A,
                                          float* buff_tau, 
                                          float* buff_beta, 
                                          float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          float* buff_u, fla_dim_t inc_u, 
                                          float* buff_a, fla_dim_t inc_a, 
                                          float* buff_y, fla_dim_t inc_y, 
                                          float* buff_w, fla_dim_t inc_w );
FLA_Error FLA_Fused_Ahx_Axpy_Ax_opd_var1( fla_dim_t m_A,
                                          fla_dim_t n_A,
                                          double* buff_tau, 
                                          double* buff_beta, 
                                          double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          double* buff_u, fla_dim_t inc_u, 
                                          double* buff_a, fla_dim_t inc_a, 
                                          double* buff_y, fla_dim_t inc_y, 
                                          double* buff_w, fla_dim_t inc_w );
FLA_Error FLA_Fused_Ahx_Axpy_Ax_opc_var1( fla_dim_t m_A,
                                          fla_dim_t n_A,
                                          scomplex* buff_tau, 
                                          scomplex* buff_beta, 
                                          scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          scomplex* buff_u, fla_dim_t inc_u, 
                                          scomplex* buff_a, fla_dim_t inc_a, 
                                          scomplex* buff_y, fla_dim_t inc_y, 
                                          scomplex* buff_w, fla_dim_t inc_w );
FLA_Error FLA_Fused_Ahx_Axpy_Ax_opz_var1( fla_dim_t m_A,
                                          fla_dim_t n_A,
                                          dcomplex* buff_tau, 
                                          dcomplex* buff_beta, 
                                          dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                          dcomplex* buff_u, fla_dim_t inc_u, 
                                          dcomplex* buff_a, fla_dim_t inc_a, 
                                          dcomplex* buff_y, fla_dim_t inc_y, 
                                          dcomplex* buff_w, fla_dim_t inc_w );

FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_opt_var1( FLA_Obj alpha, FLA_Obj tau, FLA_Obj u, FLA_Obj y, FLA_Obj z, FLA_Obj v, FLA_Obj A, FLA_Obj up, FLA_Obj a, FLA_Obj w );
FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_ops_var1( fla_dim_t m_A,
                                                fla_dim_t n_A,
                                                float* buff_tau, 
                                                float* buff_alpha, 
                                                float* buff_u, fla_dim_t inc_u, 
                                                float* buff_y, fla_dim_t inc_y, 
                                                float* buff_z, fla_dim_t inc_z, 
                                                float* buff_v, fla_dim_t inc_v, 
                                                float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                                float* buff_up, fla_dim_t inc_up, 
                                                float* buff_a, fla_dim_t inc_a, 
                                                float* buff_w, fla_dim_t inc_w );
FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_opd_var1( fla_dim_t m_A,
                                                fla_dim_t n_A,
                                                double* buff_tau, 
                                                double* buff_alpha, 
                                                double* buff_u, fla_dim_t inc_u, 
                                                double* buff_y, fla_dim_t inc_y, 
                                                double* buff_z, fla_dim_t inc_z, 
                                                double* buff_v, fla_dim_t inc_v, 
                                                double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                                double* buff_up, fla_dim_t inc_up, 
                                                double* buff_a, fla_dim_t inc_a, 
                                                double* buff_w, fla_dim_t inc_w );
FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_opc_var1( fla_dim_t m_A,
                                                fla_dim_t n_A,
                                                scomplex* buff_tau, 
                                                scomplex* buff_alpha, 
                                                scomplex* buff_u, fla_dim_t inc_u, 
                                                scomplex* buff_y, fla_dim_t inc_y, 
                                                scomplex* buff_z, fla_dim_t inc_z, 
                                                scomplex* buff_v, fla_dim_t inc_v, 
                                                scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                                scomplex* buff_up, fla_dim_t inc_up, 
                                                scomplex* buff_a, fla_dim_t inc_a, 
                                                scomplex* buff_w, fla_dim_t inc_w );
FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_opz_var1( fla_dim_t m_A,
                                                fla_dim_t n_A,
                                                dcomplex* buff_tau, 
                                                dcomplex* buff_alpha, 
                                                dcomplex* buff_u, fla_dim_t inc_u, 
                                                dcomplex* buff_y, fla_dim_t inc_y, 
                                                dcomplex* buff_z, fla_dim_t inc_z, 
                                                dcomplex* buff_v, fla_dim_t inc_v, 
                                                dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                                dcomplex* buff_up, fla_dim_t inc_up, 
                                                dcomplex* buff_a, fla_dim_t inc_a, 
                                                dcomplex* buff_w, fla_dim_t inc_w );

FLA_Error FLA_Fused_UYx_ZVx_opt_var1( FLA_Obj delta, FLA_Obj a, FLA_Obj U, FLA_Obj Y, FLA_Obj Z, FLA_Obj V, FLA_Obj A, FLA_Obj temp, FLA_Obj t, FLA_Obj w, FLA_Obj al );
FLA_Error FLA_Fused_UYx_ZVx_ops_var1( fla_dim_t m_U,
                                      fla_dim_t n_U,
                                      fla_dim_t m_V,
                                      fla_dim_t n_V,
                                      float* buff_delta, 
                                      float* buff_U, fla_dim_t rs_U, fla_dim_t cs_U, 
                                      float* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                      float* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                      float* buff_V, fla_dim_t rs_V, fla_dim_t cs_V, 
                                      float* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                      float* buff_temp, fla_dim_t inc_temp, 
                                      float* buff_t, fla_dim_t inc_t, 
                                      float* buff_a, fla_dim_t inc_a, 
                                      float* buff_w, fla_dim_t inc_w, 
                                      float* buff_al, fla_dim_t inc_al );
FLA_Error FLA_Fused_UYx_ZVx_opd_var1( fla_dim_t m_U,
                                      fla_dim_t n_U,
                                      fla_dim_t m_V,
                                      fla_dim_t n_V,
                                      double* buff_delta, 
                                      double* buff_U, fla_dim_t rs_U, fla_dim_t cs_U, 
                                      double* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                      double* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                      double* buff_V, fla_dim_t rs_V, fla_dim_t cs_V, 
                                      double* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                      double* buff_temp, fla_dim_t inc_temp, 
                                      double* buff_t, fla_dim_t inc_t, 
                                      double* buff_a, fla_dim_t inc_a, 
                                      double* buff_w, fla_dim_t inc_w, 
                                      double* buff_al, fla_dim_t inc_al );
FLA_Error FLA_Fused_UYx_ZVx_opc_var1( fla_dim_t m_U,
                                      fla_dim_t n_U,
                                      fla_dim_t m_V,
                                      fla_dim_t n_V,
                                      scomplex* buff_delta, 
                                      scomplex* buff_U, fla_dim_t rs_U, fla_dim_t cs_U, 
                                      scomplex* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                      scomplex* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                      scomplex* buff_V, fla_dim_t rs_V, fla_dim_t cs_V, 
                                      scomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                      scomplex* buff_temp, fla_dim_t inc_temp, 
                                      scomplex* buff_t, fla_dim_t inc_t, 
                                      scomplex* buff_a, fla_dim_t inc_a, 
                                      scomplex* buff_w, fla_dim_t inc_w, 
                                      scomplex* buff_al, fla_dim_t inc_al );
FLA_Error FLA_Fused_UYx_ZVx_opz_var1( fla_dim_t m_U,
                                      fla_dim_t n_U,
                                      fla_dim_t m_V,
                                      fla_dim_t n_V,
                                      dcomplex* buff_delta, 
                                      dcomplex* buff_U, fla_dim_t rs_U, fla_dim_t cs_U, 
                                      dcomplex* buff_Y, fla_dim_t rs_Y, fla_dim_t cs_Y, 
                                      dcomplex* buff_Z, fla_dim_t rs_Z, fla_dim_t cs_Z, 
                                      dcomplex* buff_V, fla_dim_t rs_V, fla_dim_t cs_V, 
                                      dcomplex* buff_A, fla_dim_t rs_A, fla_dim_t cs_A, 
                                      dcomplex* buff_temp, fla_dim_t inc_temp, 
                                      dcomplex* buff_t, fla_dim_t inc_t, 
                                      dcomplex* buff_a, fla_dim_t inc_a, 
                                      dcomplex* buff_w, fla_dim_t inc_w, 
                                      dcomplex* buff_al, fla_dim_t inc_al );
