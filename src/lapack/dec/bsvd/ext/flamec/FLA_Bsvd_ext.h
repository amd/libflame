/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- FLA_Bsvd_ext_opt_var1() ---------------------------------------------------

FLA_Error FLA_Bsvd_ext_opt_var1( fla_dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj H, 
                                 FLA_Svd_type jobu, FLA_Obj U, 
                                 FLA_Svd_type jobv, FLA_Obj V, 
                                 FLA_Bool apply_Uh2C, FLA_Obj C,
                                 fla_dim_t b_alg );
FLA_Error FLA_Bsvd_ext_ops_var1( fla_dim_t       m_d,
                                 fla_dim_t       m_U,
                                 fla_dim_t       m_V,
                                 fla_dim_t       m_C,
                                 fla_dim_t       n_C,
                                 fla_dim_t       n_GH,
                                 fla_dim_t       n_iter_max,
                                 float*    buff_d, fla_dim_t inc_d, 
                                 float*    buff_e, fla_dim_t inc_e,
                                 scomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                                 scomplex* buff_H, fla_dim_t rs_H, fla_dim_t cs_H,
                                 float*    buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                                 float*    buff_V, fla_dim_t rs_V, fla_dim_t cs_V,
                                 float*    buff_C, fla_dim_t rs_C, fla_dim_t cs_C,
                                 fla_dim_t       b_alg );
FLA_Error FLA_Bsvd_ext_opd_var1( fla_dim_t       m_d,
                                 fla_dim_t       m_U,
                                 fla_dim_t       m_V,
                                 fla_dim_t       m_C,
                                 fla_dim_t       n_C,
                                 fla_dim_t       n_GH,
                                 fla_dim_t       n_iter_max,
                                 double*   buff_d, fla_dim_t inc_d, 
                                 double*   buff_e, fla_dim_t inc_e,
                                 dcomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                                 dcomplex* buff_H, fla_dim_t rs_H, fla_dim_t cs_H,
                                 double*   buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                                 double*   buff_V, fla_dim_t rs_V, fla_dim_t cs_V,
                                 double*   buff_C, fla_dim_t rs_C, fla_dim_t cs_C,
                                 fla_dim_t       b_alg );
FLA_Error FLA_Bsvd_ext_opc_var1( fla_dim_t       m_d,
                                 fla_dim_t       m_U,
                                 fla_dim_t       m_V,
                                 fla_dim_t       m_C,
                                 fla_dim_t       n_C,
                                 fla_dim_t       n_GH,
                                 fla_dim_t       n_iter_max,
                                 float*    buff_d, fla_dim_t inc_d, 
                                 float*    buff_e, fla_dim_t inc_e,
                                 scomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                                 scomplex* buff_H, fla_dim_t rs_H, fla_dim_t cs_H,
                                 scomplex* buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                                 scomplex* buff_V, fla_dim_t rs_V, fla_dim_t cs_V,
                                 scomplex* buff_C, fla_dim_t rs_C, fla_dim_t cs_C,
                                 fla_dim_t       b_alg );
FLA_Error FLA_Bsvd_ext_opz_var1( fla_dim_t       m_d,
                                 fla_dim_t       m_U,
                                 fla_dim_t       m_V,
                                 fla_dim_t       m_C,
                                 fla_dim_t       n_C,
                                 fla_dim_t       n_GH,
                                 fla_dim_t       n_iter_max,
                                 double*   buff_d, fla_dim_t inc_d, 
                                 double*   buff_e, fla_dim_t inc_e,
                                 dcomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                                 dcomplex* buff_H, fla_dim_t rs_H, fla_dim_t cs_H,
                                 dcomplex* buff_U, fla_dim_t rs_U, fla_dim_t cs_U,
                                 dcomplex* buff_V, fla_dim_t rs_V, fla_dim_t cs_V,
                                 dcomplex* buff_C, fla_dim_t rs_C, fla_dim_t cs_C,
                                 fla_dim_t       b_alg );

