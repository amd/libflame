/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/* Test API function declaration */
void fla_test_steqr(integer argc, char **argv, test_params_t *params);
void fla_test_stevd(integer argc, char **argv, test_params_t *params);
void fla_test_geev(integer argc, char **argv, test_params_t *params);
void fla_test_geevx(integer argc, char **argv, test_params_t *params);
void fla_test_gesdd(integer argc, char **argv, test_params_t *params);
void fla_test_geqrf(integer argc, char **argv, test_params_t *params);
void fla_test_gerqf(integer argc, char **argv, test_params_t *params);
void fla_test_gerq2(integer argc, char **argv, test_params_t *params);
void fla_test_gelqf(integer argc, char **argv, test_params_t *params);
void fla_test_geqp3(integer argc, char **argv, test_params_t *params);
void fla_test_potrf(integer argc, char **argv, test_params_t *params);
void fla_test_getrf(integer argc, char **argv, test_params_t *params);
void fla_test_getri(integer argc, char **argv, test_params_t *params);
void fla_test_orgqr(integer argc, char **argv, test_params_t *params);
void fla_test_potrs(integer argc, char **argv, test_params_t *params);
void fla_test_getrs(integer argc, char **argv, test_params_t *params);
void fla_test_syevd(integer argc, char **argv, test_params_t *params);
void fla_test_gesvd(integer argc, char **argv, test_params_t *params);
void fla_test_ggevx(integer argc, char **argv, test_params_t *params);
void fla_test_gesv(integer argc, char **argv, test_params_t *params);
void fla_test_ggev(integer argc, char **argv, test_params_t *params);
void fla_test_stedc(integer argc, char **argv, test_params_t *params);
void fla_test_hseqr(integer argc, char **argv, test_params_t *params);
void fla_test_syev(integer argc, char **argv, test_params_t *params);
void fla_test_spffrt2(integer argc, char **argv, test_params_t *params);
void fla_test_spffrtx(integer argc, char **argv, test_params_t *params);
void fla_test_gehrd(integer argc, char **argv, test_params_t *params);
void fla_test_hgeqz(integer argc, char **argv, test_params_t *params);
void fla_test_gghrd(integer argc, char **argv, test_params_t *params);
void fla_test_rot(integer argc, char **argv, test_params_t *params);
void fla_test_lartg(integer argc, char **argv, test_params_t *params);
void fla_test_org2r(integer argc, char **argv, test_params_t *params);
void fla_test_syevx(integer argc, char **argv, test_params_t *params);
void fla_test_gels(integer argc, char **argv, test_params_t *params);
void fla_test_larfg(integer argc, char **argv, test_params_t *params);
void fla_test_gtsv(integer argc, char **argv, test_params_t *params);
void fla_test_gesvdx(integer argc, char **argv, test_params_t *params);
void fla_test_gbtrf(integer argc, char **argv, test_params_t *params);
void fla_test_gbtrs(integer argc, char **argv, test_params_t *params);
void fla_test_gelsd(integer argc, char **argv, test_params_t *params);
void fla_test_gelss(integer argc, char **argv, test_params_t *params);
void fla_test_sytrf(integer argc, char **argv, test_params_t *params);
void fla_test_hetrf_rook(integer argc, char **argv, test_params_t *params);
void fla_test_sytrf_rook(integer argc, char **argv, test_params_t *params);
void fla_test_larf(integer argc, char **argv, test_params_t *params);
void fla_test_sygvd(integer argc, char **argv, test_params_t *params);
void fla_test_lange(integer argc, char **argv, test_params_t *params);
void fla_test_gecon(integer argc, char **argv, test_params_t *params);
void fla_test_getrfnp(integer argc, char **argv, test_params_t *params);
void fla_test_getrfnpi(integer argc, char **argv, test_params_t *params);

#define LIN_ID 0
#define EIG_ID 1
#define SVD_ID 2
#define AUX_ID 3

/* Add test api function call entry below */
OPERATIONS API_test_functions[]
    = {{LIN_ID, "gbtrf", fla_test_gbtrf},     {LIN_ID, "gbtrs", fla_test_gbtrs},
       {LIN_ID, "orgqr", fla_test_orgqr},     {LIN_ID, "ungqr", fla_test_orgqr},
       {LIN_ID, "potrs", fla_test_potrs},     {EIG_ID, "geev", fla_test_geev},
       {EIG_ID, "geevx", fla_test_geevx},     {SVD_ID, "gesdd", fla_test_gesdd},
       {LIN_ID, "potrf", fla_test_potrf},     {LIN_ID, "geqrf", fla_test_geqrf},
       {LIN_ID, "gerqf", fla_test_gerqf},     {LIN_ID, "gerq2", fla_test_gerq2},
       {LIN_ID, "gelqf", fla_test_gelqf},     {LIN_ID, "geqp3", fla_test_geqp3},
       {LIN_ID, "getrf", fla_test_getrf},     {LIN_ID, "getri", fla_test_getri},
       {LIN_ID, "getrs", fla_test_getrs},     {EIG_ID, "syevd", fla_test_syevd},
       {EIG_ID, "heevd", fla_test_syevd},     {SVD_ID, "gesvd", fla_test_gesvd},
       {EIG_ID, "ggevx", fla_test_ggevx},     {LIN_ID, "gesv", fla_test_gesv},
       {EIG_ID, "ggev", fla_test_ggev},       {EIG_ID, "steqr", fla_test_steqr},
       {EIG_ID, "stevd", fla_test_stevd},     {EIG_ID, "stedc", fla_test_stedc},
       {EIG_ID, "hseqr", fla_test_hseqr},     {EIG_ID, "syev", fla_test_syev},
       {EIG_ID, "heev", fla_test_syev},       {LIN_ID, "spffrt2", fla_test_spffrt2},
       {LIN_ID, "spffrtx", fla_test_spffrtx}, {LIN_ID, "gehrd", fla_test_gehrd},
       {LIN_ID, "gghrd", fla_test_gghrd},     {EIG_ID, "hgeqz", fla_test_hgeqz},
       {AUX_ID, "rot", fla_test_rot},         {AUX_ID, "lartg", fla_test_lartg},
       {LIN_ID, "org2r", fla_test_org2r},     {LIN_ID, "ung2r", fla_test_org2r},
       {EIG_ID, "syevx", fla_test_syevx},     {EIG_ID, "heevx", fla_test_syevx},
       {LIN_ID, "gtsv", fla_test_gtsv},       {LIN_ID, "gels", fla_test_gels},
       {AUX_ID, "larfg", fla_test_larfg},     {SVD_ID, "gesvdx", fla_test_gesvdx},
       {LIN_ID, "gelsd", fla_test_gelsd},     {LIN_ID, "gelss", fla_test_gelss},
       {LIN_ID, "sytrf", fla_test_sytrf},     {AUX_ID, "larf", fla_test_larf},
       {EIG_ID, "sygvd", fla_test_sygvd},     {AUX_ID, "hegvd", fla_test_sygvd},
       {LIN_ID, "hetrf_rook", fla_test_hetrf_rook}, {LIN_ID, "sytrf_rook", fla_test_sytrf_rook},
       {AUX_ID, "lange", fla_test_lange},     {LIN_ID, "gecon", fla_test_gecon},
       {LIN_ID, "getrfnp", fla_test_getrfnp}, {LIN_ID, "getrfnpi", fla_test_getrfnpi}};

/* Add test API's group entry below */
char *API_test_group[] = {"LIN", "EIG", "SVD", "AUX"};
