/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Fused Level-1 BLAS-like prototypes --------------------------------------

// --- axmyv2 ---

void bl1_saxmyv2( conj1_t conjx, fla_dim_t n, float*    alpha, float*    beta, float*    x, fla_dim_t inc_x, float*    y, fla_dim_t inc_y, float*    z, fla_dim_t inc_z );
void bl1_daxmyv2( conj1_t conjx, fla_dim_t n, double*   alpha, double*   beta, double*   x, fla_dim_t inc_x, double*   y, fla_dim_t inc_y, double*   z, fla_dim_t inc_z );
void bl1_caxmyv2( conj1_t conjx, fla_dim_t n, scomplex* alpha, scomplex* beta, scomplex* x, fla_dim_t inc_x, scomplex* y, fla_dim_t inc_y, scomplex* z, fla_dim_t inc_z );
void bl1_zaxmyv2( conj1_t conjx, fla_dim_t n, dcomplex* alpha, dcomplex* beta, dcomplex* x, fla_dim_t inc_x, dcomplex* y, fla_dim_t inc_y, dcomplex* z, fla_dim_t inc_z );

// --- axpyv2b ---

void bl1_saxpyv2b( fla_dim_t n, float*    beta1, float*    beta2, float*    a1, fla_dim_t inc_a1, float*    a2, fla_dim_t inc_a2, float*    w, fla_dim_t inc_w );
void bl1_daxpyv2b( fla_dim_t n, double*   beta1, double*   beta2, double*   a1, fla_dim_t inc_a1, double*   a2, fla_dim_t inc_a2, double*   w, fla_dim_t inc_w );
void bl1_caxpyv2b( fla_dim_t n, scomplex* beta1, scomplex* beta2, scomplex* a1, fla_dim_t inc_a1, scomplex* a2, fla_dim_t inc_a2, scomplex* w, fla_dim_t inc_w );
void bl1_zaxpyv2b( fla_dim_t n, dcomplex* beta1, dcomplex* beta2, dcomplex* a1, fla_dim_t inc_a1, dcomplex* a2, fla_dim_t inc_a2, dcomplex* w, fla_dim_t inc_w );

// --- axpyv3b ---

void bl1_saxpyv3b( fla_dim_t n, float*    beta1, float*    beta2, float*    beta3, float*    a1, fla_dim_t inc_a1, float*    a2, fla_dim_t inc_a2, float*    a3, fla_dim_t inc_a3, float*    w, fla_dim_t inc_w );
void bl1_daxpyv3b( fla_dim_t n, double*   beta1, double*   beta2, double*   beta3, double*   a1, fla_dim_t inc_a1, double*   a2, fla_dim_t inc_a2, double*   a3, fla_dim_t inc_a3, double*   w, fla_dim_t inc_w );
void bl1_caxpyv3b( fla_dim_t n, scomplex* beta1, scomplex* beta2, scomplex* beta3, scomplex* a1, fla_dim_t inc_a1, scomplex* a2, fla_dim_t inc_a2, scomplex* a3, fla_dim_t inc_a3, scomplex* w, fla_dim_t inc_w );
void bl1_zaxpyv3b( fla_dim_t n, dcomplex* beta1, dcomplex* beta2, dcomplex* beta3, dcomplex* a1, fla_dim_t inc_a1, dcomplex* a2, fla_dim_t inc_a2, dcomplex* a3, fla_dim_t inc_a3, dcomplex* w, fla_dim_t inc_w );

// --- axpyv2bdotaxpy ---

void bl1_saxpyv2bdotaxpy( fla_dim_t n, float*    beta, float*    u, fla_dim_t inc_u, float*    gamma, float*    z, fla_dim_t inc_z, float*    a, fla_dim_t inc_a, float*    x, fla_dim_t inc_x, float*    kappa, float*    rho, float*    w, fla_dim_t inc_w );
void bl1_daxpyv2bdotaxpy( fla_dim_t n, double*   beta, double*   u, fla_dim_t inc_u, double*   gamma, double*   z, fla_dim_t inc_z, double*   a, fla_dim_t inc_a, double*   x, fla_dim_t inc_x, double*   kappa, double*   rho, double*   w, fla_dim_t inc_w );
void bl1_caxpyv2bdotaxpy( fla_dim_t n, scomplex* beta, scomplex* u, fla_dim_t inc_u, scomplex* gamma, scomplex* z, fla_dim_t inc_z, scomplex* a, fla_dim_t inc_a, scomplex* x, fla_dim_t inc_x, scomplex* kappa, scomplex* rho, scomplex* w, fla_dim_t inc_w );
void bl1_zaxpyv2bdotaxpy( fla_dim_t n, dcomplex* beta, dcomplex* u, fla_dim_t inc_u, dcomplex* gamma, dcomplex* z, fla_dim_t inc_z, dcomplex* a, fla_dim_t inc_a, dcomplex* x, fla_dim_t inc_x, dcomplex* kappa, dcomplex* rho, dcomplex* w, fla_dim_t inc_w );

// --- dotsv2 ---

void bl1_sdotsv2( conj1_t conjxy, fla_dim_t n, float*    x, fla_dim_t inc_x, float*    y, fla_dim_t inc_y, float*    z, fla_dim_t inc_z, float*    beta, float*    rho_xz, float*    rho_yz );
void bl1_ddotsv2( conj1_t conjxy, fla_dim_t n, double*   x, fla_dim_t inc_x, double*   y, fla_dim_t inc_y, double*   z, fla_dim_t inc_z, double*   beta, double*   rho_xz, double*   rho_yz );
void bl1_cdotsv2( conj1_t conjxy, fla_dim_t n, scomplex* x, fla_dim_t inc_x, scomplex* y, fla_dim_t inc_y, scomplex* z, fla_dim_t inc_z, scomplex* beta, scomplex* rho_xz, scomplex* rho_yz );
void bl1_zdotsv2( conj1_t conjxy, fla_dim_t n, dcomplex* x, fla_dim_t inc_x, dcomplex* y, fla_dim_t inc_y, dcomplex* z, fla_dim_t inc_z, dcomplex* beta, dcomplex* rho_xz, dcomplex* rho_yz );

// --- dotsv3 ---

void bl1_sdotsv3( conj1_t conjxyw, fla_dim_t n, float*    x, fla_dim_t inc_x, float*    y, fla_dim_t inc_y, float*    w, fla_dim_t inc_w, float*    z, fla_dim_t inc_z, float*    beta, float*    rho_xz, float*    rho_yz, float*    rho_wz );
void bl1_ddotsv3( conj1_t conjxyw, fla_dim_t n, double*   x, fla_dim_t inc_x, double*   y, fla_dim_t inc_y, double*   w, fla_dim_t inc_w, double*   z, fla_dim_t inc_z, double*   beta, double*   rho_xz, double*   rho_yz, double*   rho_wz );
void bl1_cdotsv3( conj1_t conjxyw, fla_dim_t n, scomplex* x, fla_dim_t inc_x, scomplex* y, fla_dim_t inc_y, scomplex* w, fla_dim_t inc_w, scomplex* z, fla_dim_t inc_z, scomplex* beta, scomplex* rho_xz, scomplex* rho_yz, scomplex* rho_wz );
void bl1_zdotsv3( conj1_t conjxyw, fla_dim_t n, dcomplex* x, fla_dim_t inc_x, dcomplex* y, fla_dim_t inc_y, dcomplex* w, fla_dim_t inc_w, dcomplex* z, fla_dim_t inc_z, dcomplex* beta, dcomplex* rho_xz, dcomplex* rho_yz, dcomplex* rho_wz );

// --- dotaxpy ---

void bl1_sdotaxpy( fla_dim_t n, float*    a, fla_dim_t inc_a, float*    x, fla_dim_t inc_x, float*    kappa, float*    rho, float*    w, fla_dim_t inc_w );
void bl1_ddotaxpy( fla_dim_t n, double*   a, fla_dim_t inc_a, double*   x, fla_dim_t inc_x, double*   kappa, double*   rho, double*   w, fla_dim_t inc_w );
void bl1_cdotaxpy( fla_dim_t n, scomplex* a, fla_dim_t inc_a, scomplex* x, fla_dim_t inc_x, scomplex* kappa, scomplex* rho, scomplex* w, fla_dim_t inc_w );
void bl1_zdotaxpy( fla_dim_t n, dcomplex* a, fla_dim_t inc_a, dcomplex* x, fla_dim_t inc_x, dcomplex* kappa, dcomplex* rho, dcomplex* w, fla_dim_t inc_w );

// --- dotaxmyv2 ---

void bl1_sdotaxmyv2( fla_dim_t n, float*    alpha, float*    beta, float*    x, fla_dim_t inc_x, float*    u, fla_dim_t inc_u, float*    rho, float*    y, fla_dim_t inc_y, float*    z, fla_dim_t inc_z );
void bl1_ddotaxmyv2( fla_dim_t n, double*   alpha, double*   beta, double*   x, fla_dim_t inc_x, double*   u, fla_dim_t inc_u, double*   rho, double*   y, fla_dim_t inc_y, double*   z, fla_dim_t inc_z );
void bl1_cdotaxmyv2( fla_dim_t n, scomplex* alpha, scomplex* beta, scomplex* x, fla_dim_t inc_x, scomplex* u, fla_dim_t inc_u, scomplex* rho, scomplex* y, fla_dim_t inc_y, scomplex* z, fla_dim_t inc_z );
void bl1_zdotaxmyv2( fla_dim_t n, dcomplex* alpha, dcomplex* beta, dcomplex* x, fla_dim_t inc_x, dcomplex* u, fla_dim_t inc_u, dcomplex* rho, dcomplex* y, fla_dim_t inc_y, dcomplex* z, fla_dim_t inc_z );

// --- dotv2axpyv2b ---

void bl1_sdotv2axpyv2b( fla_dim_t n, float*    a1, fla_dim_t inc_a1, float*    a2, fla_dim_t inc_a2, float*    x,  fla_dim_t inc_x, float*    kappa1, float*    kappa2, float*    rho1, float*    rho2, float*    w, fla_dim_t inc_w );
void bl1_ddotv2axpyv2b( fla_dim_t n, double*   a1, fla_dim_t inc_a1, double*   a2, fla_dim_t inc_a2, double*   x,  fla_dim_t inc_x, double*   kappa1, double*   kappa2, double*   rho1, double*   rho2, double*   w, fla_dim_t inc_w );
void bl1_cdotv2axpyv2b( fla_dim_t n, scomplex* a1, fla_dim_t inc_a1, scomplex* a2, fla_dim_t inc_a2, scomplex* x,  fla_dim_t inc_x, scomplex* kappa1, scomplex* kappa2, scomplex* rho1, scomplex* rho2, scomplex* w, fla_dim_t inc_w );
void bl1_zdotv2axpyv2b( fla_dim_t n, dcomplex* a1, fla_dim_t inc_a1, dcomplex* a2, fla_dim_t inc_a2, dcomplex* x,  fla_dim_t inc_x, dcomplex* kappa1, dcomplex* kappa2, dcomplex* rho1, dcomplex* rho2, dcomplex* w, fla_dim_t inc_w );

// --- axpyv2bdots ---

void bl1_zaxpyv2bdots( fla_dim_t       n,
                       dcomplex* alpha1,
                       dcomplex* alpha2,
                       dcomplex* x1, fla_dim_t inc_x1,
                       dcomplex* x2, fla_dim_t inc_x2,
                       dcomplex* y,  fla_dim_t inc_y,
                       dcomplex* u,  fla_dim_t inc_u,
                       dcomplex* beta,
                       dcomplex* rho );
