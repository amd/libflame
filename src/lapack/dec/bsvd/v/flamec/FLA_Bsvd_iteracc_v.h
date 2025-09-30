/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- FLA_Bsvd_iteracc_v_opt_var1() -------------------------------------------

FLA_Error FLA_Bsvd_iteracc_v_ops_var1( fla_dim_t       m_A,
                                       fla_dim_t       n_GH,
                                       fla_dim_t       ijTL,
                                       float     tol,
                                       float     thresh,
                                       float*    buff_d, fla_dim_t inc_d, 
                                       float*    buff_e, fla_dim_t inc_e,
                                       scomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                                       scomplex* buff_H, fla_dim_t rs_H, fla_dim_t cs_H,
                                       fla_dim_t*      n_iter_perf );
FLA_Error FLA_Bsvd_iteracc_v_opd_var1( fla_dim_t       m_A,
                                       fla_dim_t       n_GH,
                                       fla_dim_t       ijTL,
                                       double    tol,
                                       double    thresh,
                                       double*   buff_d, fla_dim_t inc_d, 
                                       double*   buff_e, fla_dim_t inc_e,
                                       dcomplex* buff_G, fla_dim_t rs_G, fla_dim_t cs_G,
                                       dcomplex* buff_H, fla_dim_t rs_H, fla_dim_t cs_H,
                                       fla_dim_t*      n_iter_perf );
