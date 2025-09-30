/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Name-mangling macro definitions -----------------------------------------

// --- Define Fortran name-mangling macro --------------------------
#define AOCL_LAPACK_FUNC(name,NAME) aocl_lapack_ ## name

#define F77_cgelst AOCL_LAPACK_FUNC( cgelst , CGELST )
#define F77_clatrs3 AOCL_LAPACK_FUNC( clatrs3 , CLATRS3 )
#define F77_ctrsyl3 AOCL_LAPACK_FUNC( ctrsyl3 , CTRSYL3 )
#define F77_dlarmm F77_FUNC( dlarmm , DLARMM )
#define F77_dgelst AOCL_LAPACK_FUNC( dgelst , DGELST )
#define F77_dlatrs3 AOCL_LAPACK_FUNC( dlatrs3 , DLATRS3 )
#define F77_dtrsyl3 AOCL_LAPACK_FUNC( dtrsyl3 , DTRSYL3 )
#define F77_slarmm F77_FUNC( slarmm , SLARMM )
#define F77_sgelst AOCL_LAPACK_FUNC( sgelst , SGELST )
#define F77_slatrs3 AOCL_LAPACK_FUNC( slatrs3 , SLATRS3 )
#define F77_strsyl3 AOCL_LAPACK_FUNC( strsyl3 , STRSYL3 )
#define F77_zgelst AOCL_LAPACK_FUNC( zgelst , ZGELST )
#define F77_zlatrs3 AOCL_LAPACK_FUNC( zlatrs3 , ZLATRS3 )
#define F77_ztrsyl3 AOCL_LAPACK_FUNC( ztrsyl3 , ZTRSYL3 )
#define F77_dlamch F77_FUNC( dlamch , DLAMCH )
#define F77_dlamc3 F77_FUNC( dlamc3 , DLAMC3 )
#define F77_dladiv2 F77_FUNC( dladiv2 , DLADIV2 )
#define F77_dladiv1 F77_FUNC(dladiv1 , DLADIV1 )
#define F77_sladiv2 F77_FUNC( sladiv2 , SLADIV2 )
#define F77_sladiv1 F77_FUNC( sladiv1 , SLADIV1 )
#define F77_iparmq AOCL_LAPACK_FUNC( iparmq, IPARMQ )
#define F77_ilaenv AOCL_LAPACK_FUNC( ilaenv , ILAENV )
#define F77_ieeeck AOCL_LAPACK_FUNC( ieeeck , IEEECK )
#define F77_lsamen AOCL_LAPACK_FUNC( lsamen , LSAMEN )
#define F77_slamc3 F77_FUNC( slamc3 , SLAMC3 )
#define F77_slamch F77_FUNC( slamch, SLAMCH )
#define F77_cgetsqrhrt AOCL_LAPACK_FUNC( cgetsqrhrt , CGETSQRHRT )
#define F77_claqz0 AOCL_LAPACK_FUNC( claqz0 , CLAQZ0 )
#define F77_claqz1 AOCL_LAPACK_FUNC( claqz1 , CLAQZ1 )
#define F77_claqz2 AOCL_LAPACK_FUNC( claqz2 , CLAQZ2 )
#define F77_claqz3 AOCL_LAPACK_FUNC( claqz3 , CLAQZ3 )
#define F77_clarfb_gett AOCL_LAPACK_FUNC( clarfb_gett , CLARFB_GETT )
#define F77_cungtsqr_row AOCL_LAPACK_FUNC( cungtsqr_row , CUNGTSQR_ROW )
#define F77_dgetsqrhrt AOCL_LAPACK_FUNC( dgetsqrhrt , DGETSQRHRT )
#define F77_dlaqz0 AOCL_LAPACK_FUNC( dlaqz0 , DLAQZ0 )
#define F77_dlaqz1 AOCL_LAPACK_FUNC( dlaqz1 , DLAQZ1 )
#define F77_dlaqz2 AOCL_LAPACK_FUNC( dlaqz2 , DLAQZ2 )
#define F77_dlaqz3 AOCL_LAPACK_FUNC( dlaqz3 , DLAQZ3 )
#define F77_dlaqz4 AOCL_LAPACK_FUNC( dlaqz4 , DLAQZ4 )
#define F77_dlarfb_gett AOCL_LAPACK_FUNC( dlarfb_gett , DLARFB_GETT )
#define F77_dorgtsqr_row AOCL_LAPACK_FUNC( dorgtsqr_row , DORGTSQR_ROW )
#define F77_sgetsqrhrt AOCL_LAPACK_FUNC( sgetsqrhrt , SGETSQRHRT )
#define F77_slaqz0 AOCL_LAPACK_FUNC( slaqz0 , SLAQZ0 )
#define F77_slaqz1 AOCL_LAPACK_FUNC( slaqz1 , SLAQZ1 )
#define F77_slaqz2 AOCL_LAPACK_FUNC( slaqz2 , SLAQZ2 )
#define F77_slaqz3 AOCL_LAPACK_FUNC( slaqz3 , SLAQZ3 )
#define F77_slaqz4 AOCL_LAPACK_FUNC( slaqz4 , SLAQZ4 )
#define F77_slarfb_gett AOCL_LAPACK_FUNC( slarfb_gett , SLARFB_GETT )
#define F77_sorgtsqr_row AOCL_LAPACK_FUNC( sorgtsqr_row , SORGTSQR_ROW )
#define F77_zgetsqrhrt AOCL_LAPACK_FUNC( zgetsqrhrt , ZGETSQRHRT )
#define F77_zlaqz0 AOCL_LAPACK_FUNC( zlaqz0 , ZLAQZ0 )
#define F77_zlaqz1 AOCL_LAPACK_FUNC( zlaqz1 , ZLAQZ1 )
#define F77_zlaqz2 AOCL_LAPACK_FUNC( zlaqz2 , ZLAQZ2 )
#define F77_zlaqz3 AOCL_LAPACK_FUNC( zlaqz3 , ZLAQZ3 )
#define F77_zlarfb_gett AOCL_LAPACK_FUNC( zlarfb_gett , ZLARFB_GETT )
#define F77_zungtsqr_row AOCL_LAPACK_FUNC( zungtsqr_row , ZUNGTSQR_ROW )
#define F77_cbbcsd AOCL_LAPACK_FUNC( cbbcsd , CBBCSD )
#define F77_cbdsqr AOCL_LAPACK_FUNC( cbdsqr , CBDSQR )
#define F77_cgbbrd AOCL_LAPACK_FUNC( cgbbrd , CGBBRD )
#define F77_cgbcon AOCL_LAPACK_FUNC( cgbcon , CGBCON )
#define F77_cgbequb AOCL_LAPACK_FUNC( cgbequb , CGBEQUB )
#define F77_cgbequ AOCL_LAPACK_FUNC( cgbequ , CGBEQU )
#define F77_cgbrfs AOCL_LAPACK_FUNC( cgbrfs , CGBRFS )
#define F77_cgbrfsx F77_FUNC( cgbrfsx , CGBRFSX )
#define F77_cgbsv AOCL_LAPACK_FUNC( cgbsv , CGBSV )
#define F77_cgbsvx AOCL_LAPACK_FUNC( cgbsvx , CGBSVX )
#define F77_cgbsvxx F77_FUNC( cgbsvxx , CGBSVXX )
#define F77_cgbtf2 AOCL_LAPACK_FUNC( cgbtf2 , CGBTF2 )
#define F77_cgbtrf AOCL_LAPACK_FUNC( cgbtrf , CGBTRF )
#define F77_cgbtrs AOCL_LAPACK_FUNC( cgbtrs , CGBTRS )
#define F77_cgebak AOCL_LAPACK_FUNC( cgebak , CGEBAK )
#define F77_cgebal AOCL_LAPACK_FUNC( cgebal , CGEBAL )
#define F77_cgebd2 AOCL_LAPACK_FUNC( cgebd2 , CGEBD2 )
#define F77_cgebrd AOCL_LAPACK_FUNC( cgebrd , CGEBRD )
#define F77_cgecon AOCL_LAPACK_FUNC( cgecon , CGECON )
#define F77_cgeequb AOCL_LAPACK_FUNC( cgeequb , CGEEQUB )
#define F77_cgeequ AOCL_LAPACK_FUNC( cgeequ , CGEEQU )
#define F77_cgees AOCL_LAPACK_FUNC( cgees , CGEES )
#define F77_cgeesx AOCL_LAPACK_FUNC( cgeesx , CGEESX )
#define F77_cgeev AOCL_LAPACK_FUNC( cgeev , CGEEV )
#define F77_cgeevx AOCL_LAPACK_FUNC( cgeevx , CGEEVX )
#define F77_cgegs AOCL_LAPACK_FUNC( cgegs , CGEGS )
#define F77_cgegv AOCL_LAPACK_FUNC( cgegv , CGEGV )
#define F77_cgehd2 AOCL_LAPACK_FUNC( cgehd2 , CGEHD2 )
#define F77_cgehrd AOCL_LAPACK_FUNC( cgehrd , CGEHRD )
#define F77_cgejsv AOCL_LAPACK_FUNC( cgejsv , CGEJSV )
#define F77_cgelq2 AOCL_LAPACK_FUNC( cgelq2 , CGELQ2 )
#define F77_cgelq AOCL_LAPACK_FUNC( cgelq , CGELQ )
#define F77_cgelqf AOCL_LAPACK_FUNC( cgelqf , CGELQF )
#define F77_cgelqt3 AOCL_LAPACK_FUNC( cgelqt3 , CGELQT3 )
#define F77_cgelqt AOCL_LAPACK_FUNC( cgelqt , CGELQT )
#define F77_cgelsd AOCL_LAPACK_FUNC( cgelsd , CGELSD )
#define F77_cgels AOCL_LAPACK_FUNC( cgels , CGELS )
#define F77_cgelss AOCL_LAPACK_FUNC( cgelss , CGELSS )
#define F77_cgelsx AOCL_LAPACK_FUNC( cgelsx , CGELSX )
#define F77_cgelsy AOCL_LAPACK_FUNC( cgelsy , CGELSY )
#define F77_cgemlq AOCL_LAPACK_FUNC( cgemlq , CGEMLQ )
#define F77_cgemlqt AOCL_LAPACK_FUNC( cgemlqt , CGEMLQT )
#define F77_cgemqr AOCL_LAPACK_FUNC( cgemqr , CGEMQR )
#define F77_cgemqrt AOCL_LAPACK_FUNC( cgemqrt , CGEMQRT )
#define F77_cgeql2 AOCL_LAPACK_FUNC( cgeql2 , CGEQL2 )
#define F77_cgeqlf AOCL_LAPACK_FUNC( cgeqlf , CGEQLF )
#define F77_cgeqp3 AOCL_LAPACK_FUNC( cgeqp3 , CGEQP3 )
#define F77_cgeqpf AOCL_LAPACK_FUNC( cgeqpf , CGEQPF )
#define F77_cgeqr2 AOCL_LAPACK_FUNC( cgeqr2 , CGEQR2 )
#define F77_cgeqr2p AOCL_LAPACK_FUNC( cgeqr2p , CGEQR2P )
#define F77_cgeqr AOCL_LAPACK_FUNC( cgeqr , CGEQR )
#define F77_cgeqrf AOCL_LAPACK_FUNC( cgeqrf , CGEQRF )
#define F77_cgeqrfp AOCL_LAPACK_FUNC( cgeqrfp , CGEQRFP )
#define F77_cgeqrt2 AOCL_LAPACK_FUNC( cgeqrt2 , CGEQRT2 )
#define F77_cgeqrt3 AOCL_LAPACK_FUNC( cgeqrt3 , CGEQRT3 )
#define F77_cgeqrt AOCL_LAPACK_FUNC( cgeqrt , CGEQRT )
#define F77_cgerfs AOCL_LAPACK_FUNC( cgerfs , CGERFS )
#define F77_cgerfsx F77_FUNC( cgerfsx , CGERFSX )
#define F77_cgerq2 AOCL_LAPACK_FUNC( cgerq2 , CGERQ2 )
#define F77_cgerqf AOCL_LAPACK_FUNC( cgerqf , CGERQF )
#define F77_cgesc2 AOCL_LAPACK_FUNC( cgesc2 , CGESC2 )
#define F77_cgesdd AOCL_LAPACK_FUNC( cgesdd , CGESDD )
#define F77_cgesvd AOCL_LAPACK_FUNC( cgesvd , CGESVD )
#define F77_cgesvdq AOCL_LAPACK_FUNC( cgesvdq , CGESVDQ )
#define F77_cgesvdx AOCL_LAPACK_FUNC( cgesvdx , CGESVDX )
#define F77_cgesv AOCL_LAPACK_FUNC( cgesv , CGESV )
#define F77_cgesvj AOCL_LAPACK_FUNC( cgesvj , CGESVJ )
#define F77_cgesvx AOCL_LAPACK_FUNC( cgesvx , CGESVX )
#define F77_cgesvxx F77_FUNC( cgesvxx , CGESVXX )
#define F77_cgetc2 AOCL_LAPACK_FUNC( cgetc2 , CGETC2 )
#define F77_cgetf2 AOCL_LAPACK_FUNC( cgetf2 , CGETF2 )
#define F77_cgetrf2 AOCL_LAPACK_FUNC( cgetrf2 , CGETRF2 )
#define F77_cgetrf AOCL_LAPACK_FUNC( cgetrf , CGETRF )
#define F77_cgetri AOCL_LAPACK_FUNC( cgetri , CGETRI )
#define F77_cgetrs AOCL_LAPACK_FUNC( cgetrs , CGETRS )
#define F77_cgetsls AOCL_LAPACK_FUNC( cgetsls , CGETSLS )
#define F77_cggbak AOCL_LAPACK_FUNC( cggbak , CGGBAK )
#define F77_cggbal AOCL_LAPACK_FUNC( cggbal , CGGBAL )
#define F77_cgges3 AOCL_LAPACK_FUNC( cgges3 , CGGES3 )
#define F77_cgges AOCL_LAPACK_FUNC( cgges , CGGES )
#define F77_cggesx AOCL_LAPACK_FUNC( cggesx , CGGESX )
#define F77_cggev3 AOCL_LAPACK_FUNC( cggev3 , CGGEV3 )
#define F77_cggev AOCL_LAPACK_FUNC( cggev , CGGEV )
#define F77_cggevx AOCL_LAPACK_FUNC( cggevx , CGGEVX )
#define F77_cggglm AOCL_LAPACK_FUNC( cggglm , CGGGLM )
#define F77_cgghd3 AOCL_LAPACK_FUNC( cgghd3 , CGGHD3 )
#define F77_cgghrd AOCL_LAPACK_FUNC( cgghrd , CGGHRD )
#define F77_cgglse AOCL_LAPACK_FUNC( cgglse , CGGLSE )
#define F77_cggqrf AOCL_LAPACK_FUNC( cggqrf , CGGQRF )
#define F77_cggrqf AOCL_LAPACK_FUNC( cggrqf , CGGRQF )
#define F77_cggsvd3 AOCL_LAPACK_FUNC( cggsvd3 , CGGSVD3 )
#define F77_cggsvd AOCL_LAPACK_FUNC( cggsvd , CGGSVD )
#define F77_cggsvp3 AOCL_LAPACK_FUNC( cggsvp3 , CGGSVP3 )
#define F77_cggsvp AOCL_LAPACK_FUNC( cggsvp , CGGSVP )
#define F77_cgsvj0 AOCL_LAPACK_FUNC( cgsvj0 , CGSVJ0 )
#define F77_cgsvj1 AOCL_LAPACK_FUNC( cgsvj1 , CGSVJ1 )
#define F77_cgtcon AOCL_LAPACK_FUNC( cgtcon , CGTCON )
#define F77_cgtrfs AOCL_LAPACK_FUNC( cgtrfs , CGTRFS )
#define F77_cgtsv AOCL_LAPACK_FUNC( cgtsv , CGTSV )
#define F77_cgtsvx AOCL_LAPACK_FUNC( cgtsvx , CGTSVX )
#define F77_cgttrf AOCL_LAPACK_FUNC( cgttrf , CGTTRF )
#define F77_cgttrs AOCL_LAPACK_FUNC( cgttrs , CGTTRS )
#define F77_cgtts2 AOCL_LAPACK_FUNC( cgtts2 , CGTTS2 )
#define F77_chb2st_kernels F77_FUNC( chb2st_kernels , CHB2ST_KERNELS )
#define F77_chbev_2stage AOCL_LAPACK_FUNC( chbev_2stage , CHBEV_2STAGE )
#define F77_chbevd_2stage AOCL_LAPACK_FUNC( chbevd_2stage , CHBEVD_2STAGE )
#define F77_chbevd AOCL_LAPACK_FUNC( chbevd , CHBEVD )
#define F77_chbev AOCL_LAPACK_FUNC( chbev , CHBEV )
#define F77_chbevx_2stage AOCL_LAPACK_FUNC( chbevx_2stage , CHBEVX_2STAGE )
#define F77_chbevx AOCL_LAPACK_FUNC( chbevx , CHBEVX )
#define F77_chbgst AOCL_LAPACK_FUNC( chbgst , CHBGST )
#define F77_chbgvd AOCL_LAPACK_FUNC( chbgvd , CHBGVD )
#define F77_chbgv AOCL_LAPACK_FUNC( chbgv , CHBGV )
#define F77_chbgvx AOCL_LAPACK_FUNC( chbgvx , CHBGVX )
#define F77_chbtrd AOCL_LAPACK_FUNC( chbtrd , CHBTRD )
#define F77_checon_3 AOCL_LAPACK_FUNC( checon_3 , CHECON_3 )
#define F77_checon AOCL_LAPACK_FUNC( checon , CHECON )
#define F77_checon_rook AOCL_LAPACK_FUNC( checon_rook , CHECON_ROOK )
#define F77_cheequb AOCL_LAPACK_FUNC( cheequb , CHEEQUB )
#define F77_cheev_2stage AOCL_LAPACK_FUNC( cheev_2stage , CHEEV_2STAGE )
#define F77_cheevd_2stage AOCL_LAPACK_FUNC( cheevd_2stage , CHEEVD_2STAGE )
#define F77_cheevd AOCL_LAPACK_FUNC( cheevd , CHEEVD )
#define F77_cheev  AOCL_LAPACK_FUNC( cheev  , CHEEV  )
#define F77_cheev AOCL_LAPACK_FUNC( cheev , CHEEV )
#define F77_cheevr_2stage AOCL_LAPACK_FUNC( cheevr_2stage , CHEEVR_2STAGE )
#define F77_cheevr AOCL_LAPACK_FUNC( cheevr , CHEEVR )
#define F77_cheevx_2stage AOCL_LAPACK_FUNC( cheevx_2stage , CHEEVX_2STAGE )
#define F77_cheevx AOCL_LAPACK_FUNC( cheevx , CHEEVX )
#define F77_chegs2 AOCL_LAPACK_FUNC( chegs2 , CHEGS2 )
#define F77_chegst AOCL_LAPACK_FUNC( chegst , CHEGST )
#define F77_chegv_2stage AOCL_LAPACK_FUNC( chegv_2stage , CHEGV_2STAGE )
#define F77_chegvd AOCL_LAPACK_FUNC( chegvd , CHEGVD )
#define F77_chegv AOCL_LAPACK_FUNC( chegv , CHEGV )
#define F77_chegvx AOCL_LAPACK_FUNC( chegvx , CHEGVX )
#define F77_cherfs AOCL_LAPACK_FUNC( cherfs , CHERFS )
#define F77_cherfsx F77_FUNC( cherfsx , CHERFSX )
#define F77_chesv_aa_2stage AOCL_LAPACK_FUNC( chesv_aa_2stage , CHESV_AA_2STAGE )
#define F77_chesv_aa AOCL_LAPACK_FUNC( chesv_aa , CHESV_AA )
#define F77_chesv AOCL_LAPACK_FUNC( chesv , CHESV )
#define F77_chesv_rk AOCL_LAPACK_FUNC( chesv_rk , CHESV_RK )
#define F77_chesv_rook AOCL_LAPACK_FUNC( chesv_rook , CHESV_ROOK )
#define F77_chesvx AOCL_LAPACK_FUNC( chesvx , CHESVX )
#define F77_chesvxx F77_FUNC( chesvxx , CHESVXX )
#define F77_cheswapr AOCL_LAPACK_FUNC( cheswapr , CHESWAPR )
#define F77_chetd2 AOCL_LAPACK_FUNC( chetd2 , CHETD2 )
#define F77_chetf2 AOCL_LAPACK_FUNC( chetf2 , CHETF2 )
#define F77_chetf2_rk AOCL_LAPACK_FUNC( chetf2_rk , CHETF2_RK )
#define F77_chetf2_rook AOCL_LAPACK_FUNC( chetf2_rook , CHETF2_ROOK )
#define F77_chetrd_2stage AOCL_LAPACK_FUNC( chetrd_2stage , CHETRD_2STAGE )
#define F77_chetrd AOCL_LAPACK_FUNC( chetrd , CHETRD )
#define F77_chetrd_hb2st AOCL_LAPACK_FUNC( chetrd_hb2st , CHETRD_HB2ST )
#define F77_chetrd_he2hb AOCL_LAPACK_FUNC( chetrd_he2hb , CHETRD_HE2HB )
#define F77_chetrf_aa_2stage AOCL_LAPACK_FUNC( chetrf_aa_2stage , CHETRF_AA_2STAGE )
#define F77_chetrf_aa AOCL_LAPACK_FUNC( chetrf_aa , CHETRF_AA )
#define F77_chetrf AOCL_LAPACK_FUNC( chetrf , CHETRF )
#define F77_chetrf_rk AOCL_LAPACK_FUNC( chetrf_rk , CHETRF_RK )
#define F77_chetrf_rook AOCL_LAPACK_FUNC( chetrf_rook , CHETRF_ROOK )
#define F77_chetri2 AOCL_LAPACK_FUNC( chetri2 , CHETRI2 )
#define F77_chetri2x AOCL_LAPACK_FUNC( chetri2x , CHETRI2X )
#define F77_chetri_3 AOCL_LAPACK_FUNC( chetri_3 , CHETRI_3 )
#define F77_chetri_3x AOCL_LAPACK_FUNC( chetri_3x , CHETRI_3X )
#define F77_chetri AOCL_LAPACK_FUNC( chetri , CHETRI )
#define F77_chetri_rook AOCL_LAPACK_FUNC( chetri_rook , CHETRI_ROOK )
#define F77_chetrs2 AOCL_LAPACK_FUNC( chetrs2 , CHETRS2 )
#define F77_chetrs_3 AOCL_LAPACK_FUNC( chetrs_3 , CHETRS_3 )
#define F77_chetrs_aa_2stage AOCL_LAPACK_FUNC( chetrs_aa_2stage , CHETRS_AA_2STAGE )
#define F77_chetrs_aa AOCL_LAPACK_FUNC( chetrs_aa , CHETRS_AA )
#define F77_chetrs AOCL_LAPACK_FUNC( chetrs , CHETRS )
#define F77_chetrs_rook AOCL_LAPACK_FUNC( chetrs_rook , CHETRS_ROOK )
#define F77_chfrk AOCL_LAPACK_FUNC( chfrk , CHFRK )
#define F77_chgeqz AOCL_LAPACK_FUNC( chgeqz , CHGEQZ )
#define F77_chla_transtype F77_FUNC( chla_transtype , CHLA_TRANSTYPE )
#define F77_chpcon AOCL_LAPACK_FUNC( chpcon , CHPCON )
#define F77_chpevd AOCL_LAPACK_FUNC( chpevd , CHPEVD )
#define F77_chpev AOCL_LAPACK_FUNC( chpev , CHPEV )
#define F77_chpevx AOCL_LAPACK_FUNC( chpevx , CHPEVX )
#define F77_chpgst AOCL_LAPACK_FUNC( chpgst , CHPGST )
#define F77_chpgvd AOCL_LAPACK_FUNC( chpgvd , CHPGVD )
#define F77_chpgv AOCL_LAPACK_FUNC( chpgv , CHPGV )
#define F77_chpgvx AOCL_LAPACK_FUNC( chpgvx , CHPGVX )
#define F77_chprfs AOCL_LAPACK_FUNC( chprfs , CHPRFS )
#define F77_chpsv AOCL_LAPACK_FUNC( chpsv , CHPSV )
#define F77_chpsvx AOCL_LAPACK_FUNC( chpsvx , CHPSVX )
#define F77_chptrd AOCL_LAPACK_FUNC( chptrd , CHPTRD )
#define F77_chptrf AOCL_LAPACK_FUNC( chptrf , CHPTRF )
#define F77_chptri AOCL_LAPACK_FUNC( chptri , CHPTRI )
#define F77_chptrs AOCL_LAPACK_FUNC( chptrs , CHPTRS )
#define F77_chsein AOCL_LAPACK_FUNC( chsein , CHSEIN )
#define F77_chseqr AOCL_LAPACK_FUNC( chseqr , CHSEQR )
#define F77_clabrd AOCL_LAPACK_FUNC( clabrd , CLABRD )
#define F77_clacgv AOCL_LAPACK_FUNC( clacgv , CLACGV )
#define F77_clacn2 AOCL_LAPACK_FUNC( clacn2 , CLACN2 )
#define F77_clacon AOCL_LAPACK_FUNC( clacon , CLACON )
#define F77_clacp2 AOCL_LAPACK_FUNC( clacp2 , CLACP2 )
#define F77_clacpy AOCL_LAPACK_FUNC( clacpy , CLACPY )
#define F77_clacrm AOCL_LAPACK_FUNC( clacrm , CLACRM )
#define F77_clacrt AOCL_LAPACK_FUNC( clacrt , CLACRT )
#define F77_cladiv F77_FUNC( cladiv , CLADIV )
#define F77_claed0 AOCL_LAPACK_FUNC( claed0 , CLAED0 )
#define F77_claed7 AOCL_LAPACK_FUNC( claed7 , CLAED7 )
#define F77_claed8 AOCL_LAPACK_FUNC( claed8 , CLAED8 )
#define F77_claein AOCL_LAPACK_FUNC( claein , CLAEIN )
#define F77_claesy F77_FUNC( claesy , CLAESY )
#define F77_claev2 F77_FUNC( claev2 , CLAEV2 )
#define F77_clag2z AOCL_LAPACK_FUNC( clag2z , CLAG2Z )
#define F77_cla_gbamv F77_FUNC( cla_gbamv , CLA_GBAMV )
#define F77_cla_gbrcond_c F77_FUNC( cla_gbrcond_c , CLA_GBRCOND_C )
#define F77_cla_gbrcond_x F77_FUNC( cla_gbrcond_x , CLA_GBRCOND_X )
#define F77_cla_gbrfsx_extended F77_FUNC( cla_gbrfsx_extended , CLA_GBRFSX_EXTENDED )
#define F77_cla_gbrpvgrw F77_FUNC( cla_gbrpvgrw , CLA_GBRPVGRW )
#define F77_cla_geamv F77_FUNC( cla_geamv , CLA_GEAMV )
#define F77_cla_gercond_c F77_FUNC( cla_gercond_c , CLA_GERCOND_C )
#define F77_cla_gercond_x F77_FUNC( cla_gercond_x , CLA_GERCOND_X )
#define F77_cla_gerfsx_extended F77_FUNC( cla_gerfsx_extended , CLA_GERFSX_EXTENDED )
#define F77_cla_gerpvgrw F77_FUNC( cla_gerpvgrw , CLA_GERPVGRW )
#define F77_clags2 F77_FUNC( clags2 , CLAGS2 )
#define F77_clagtm AOCL_LAPACK_FUNC( clagtm , CLAGTM )
#define F77_cla_heamv F77_FUNC( cla_heamv , CLA_HEAMV )
#define F77_clahef_aa AOCL_LAPACK_FUNC( clahef_aa , CLAHEF_AA )
#define F77_clahef AOCL_LAPACK_FUNC( clahef , CLAHEF )
#define F77_clahef_rk AOCL_LAPACK_FUNC( clahef_rk , CLAHEF_RK )
#define F77_clahef_rook AOCL_LAPACK_FUNC( clahef_rook , CLAHEF_ROOK )
#define F77_cla_hercond_c F77_FUNC( cla_hercond_c , CLA_HERCOND_C )
#define F77_cla_hercond_x F77_FUNC( cla_hercond_x , CLA_HERCOND_X )
#define F77_cla_herfsx_extended F77_FUNC( cla_herfsx_extended , CLA_HERFSX_EXTENDED )
#define F77_cla_herpvgrw F77_FUNC( cla_herpvgrw , CLA_HERPVGRW )
#define F77_clahqr AOCL_LAPACK_FUNC( clahqr , CLAHQR )
#define F77_clahr2 AOCL_LAPACK_FUNC( clahr2 , CLAHR2 )
#define F77_clahrd AOCL_LAPACK_FUNC( clahrd , CLAHRD )
#define F77_claic1 AOCL_LAPACK_FUNC( claic1 , CLAIC1 )
#define F77_cla_lin_berr F77_FUNC( cla_lin_berr , CLA_LIN_BERR )
#define F77_clals0 AOCL_LAPACK_FUNC( clals0 , CLALS0 )
#define F77_clalsa AOCL_LAPACK_FUNC( clalsa , CLALSA )
#define F77_clalsd AOCL_LAPACK_FUNC( clalsd , CLALSD )
#define F77_clamswlq AOCL_LAPACK_FUNC( clamswlq , CLAMSWLQ )
#define F77_clamtsqr AOCL_LAPACK_FUNC( clamtsqr , CLAMTSQR )
#define F77_clangb AOCL_LAPACK_FUNC( clangb , CLANGB )
#define F77_clange AOCL_LAPACK_FUNC( clange , CLANGE )
#define F77_clangt AOCL_LAPACK_FUNC( clangt , CLANGT )
#define F77_clanhb AOCL_LAPACK_FUNC( clanhb , CLANHB )
#define F77_clanhe AOCL_LAPACK_FUNC( clanhe , CLANHE )
#define F77_clanhf AOCL_LAPACK_FUNC( clanhf , CLANHF )
#define F77_clanhp AOCL_LAPACK_FUNC( clanhp , CLANHP )
#define F77_clanhs AOCL_LAPACK_FUNC( clanhs , CLANHS )
#define F77_clanht AOCL_LAPACK_FUNC( clanht , CLANHT )
#define F77_clansb AOCL_LAPACK_FUNC( clansb , CLANSB )
#define F77_clansp AOCL_LAPACK_FUNC( clansp , CLANSP )
#define F77_clansy AOCL_LAPACK_FUNC( clansy , CLANSY )
#define F77_clantb AOCL_LAPACK_FUNC( clantb , CLANTB )
#define F77_clantp AOCL_LAPACK_FUNC( clantp , CLANTP )
#define F77_clantr AOCL_LAPACK_FUNC( clantr , CLANTR )
#define F77_clapll AOCL_LAPACK_FUNC( clapll , CLAPLL )
#define F77_clapmr AOCL_LAPACK_FUNC( clapmr , CLAPMR )
#define F77_clapmt AOCL_LAPACK_FUNC( clapmt , CLAPMT )
#define F77_cla_porcond_c F77_FUNC( cla_porcond_c , CLA_PORCOND_C )
#define F77_cla_porcond_x F77_FUNC( cla_porcond_x , CLA_PORCOND_X )
#define F77_cla_porfsx_extended F77_FUNC( cla_porfsx_extended , CLA_PORFSX_EXTENDED )
#define F77_cla_porpvgrw F77_FUNC( cla_porpvgrw , CLA_PORPVGRW )
#define F77_claqgb AOCL_LAPACK_FUNC( claqgb , CLAQGB )
#define F77_claqge AOCL_LAPACK_FUNC( claqge , CLAQGE )
#define F77_claqhb AOCL_LAPACK_FUNC( claqhb , CLAQHB )
#define F77_claqhe AOCL_LAPACK_FUNC( claqhe , CLAQHE )
#define F77_claqhp AOCL_LAPACK_FUNC( claqhp , CLAQHP )
#define F77_claqp2 AOCL_LAPACK_FUNC( claqp2 , CLAQP2 )
#define F77_claqps AOCL_LAPACK_FUNC( claqps , CLAQPS )
#define F77_claqr0 AOCL_LAPACK_FUNC( claqr0 , CLAQR0 )
#define F77_claqr1 AOCL_LAPACK_FUNC( claqr1 , CLAQR1 )
#define F77_claqr2 AOCL_LAPACK_FUNC( claqr2 , CLAQR2 )
#define F77_claqr3 AOCL_LAPACK_FUNC( claqr3 , CLAQR3 )
#define F77_claqr4 AOCL_LAPACK_FUNC( claqr4 , CLAQR4 )
#define F77_claqr5 AOCL_LAPACK_FUNC( claqr5 , CLAQR5 )
#define F77_claqsb AOCL_LAPACK_FUNC( claqsb , CLAQSB )
#define F77_claqsp AOCL_LAPACK_FUNC( claqsp , CLAQSP )
#define F77_claqsy AOCL_LAPACK_FUNC( claqsy , CLAQSY )
#define F77_clar1v AOCL_LAPACK_FUNC( clar1v , CLAR1V )
#define F77_clar2v AOCL_LAPACK_FUNC( clar2v , CLAR2V )
#define F77_clarcm AOCL_LAPACK_FUNC( clarcm , CLARCM )
#define F77_clarfb AOCL_LAPACK_FUNC( clarfb , CLARFB )
#define F77_clarf AOCL_LAPACK_FUNC( clarf , CLARF )
#define F77_clarfg AOCL_LAPACK_FUNC( clarfg , CLARFG )
#define F77_clarfgp AOCL_LAPACK_FUNC( clarfgp , CLARFGP )
#define F77_clarft AOCL_LAPACK_FUNC( clarft , CLARFT )
#define F77_clarfx AOCL_LAPACK_FUNC( clarfx , CLARFX )
#define F77_clarfy AOCL_LAPACK_FUNC( clarfy , CLARFY )
#define F77_clargv AOCL_LAPACK_FUNC( clargv , CLARGV )
#define F77_clarnv AOCL_LAPACK_FUNC( clarnv , CLARNV )
#define F77_clarrv AOCL_LAPACK_FUNC( clarrv , CLARRV )
#define F77_clarscl2 AOCL_LAPACK_FUNC( clarscl2 , CLARSCL2 )
#define F77_clartg F77_FUNC( clartg , CLARTG )
#define F77_clartv AOCL_LAPACK_FUNC( clartv , CLARTV )
#define F77_clarzb AOCL_LAPACK_FUNC( clarzb , CLARZB )
#define F77_clarz AOCL_LAPACK_FUNC( clarz , CLARZ )
#define F77_clarzt AOCL_LAPACK_FUNC( clarzt , CLARZT )
#define F77_clascl2 AOCL_LAPACK_FUNC( clascl2 , CLASCL2 )
#define F77_clascl AOCL_LAPACK_FUNC( clascl , CLASCL )
#define F77_claset AOCL_LAPACK_FUNC( claset , CLASET )
#define F77_clasr AOCL_LAPACK_FUNC( clasr , CLASR )
#define F77_classq AOCL_LAPACK_FUNC( classq , CLASSQ )
#define F77_claswlq AOCL_LAPACK_FUNC( claswlq , CLASWLQ )
#define F77_claswp AOCL_LAPACK_FUNC( claswp , CLASWP )
#define F77_cla_syamv F77_FUNC( cla_syamv , CLA_SYAMV )
#define F77_clasyf_aa AOCL_LAPACK_FUNC( clasyf_aa , CLASYF_AA )
#define F77_clasyf AOCL_LAPACK_FUNC( clasyf , CLASYF )
#define F77_clasyf_rk AOCL_LAPACK_FUNC( clasyf_rk , CLASYF_RK )
#define F77_clasyf_rook AOCL_LAPACK_FUNC( clasyf_rook , CLASYF_ROOK )
#define F77_cla_syrcond_c F77_FUNC( cla_syrcond_c , CLA_SYRCOND_C )
#define F77_cla_syrcond_x F77_FUNC( cla_syrcond_x , CLA_SYRCOND_X )
#define F77_cla_syrfsx_extended F77_FUNC( cla_syrfsx_extended , CLA_SYRFSX_EXTENDED )
#define F77_cla_syrpvgrw F77_FUNC( cla_syrpvgrw , CLA_SYRPVGRW )
#define F77_clatbs AOCL_LAPACK_FUNC( clatbs , CLATBS )
#define F77_clatdf AOCL_LAPACK_FUNC( clatdf , CLATDF )
#define F77_clatps AOCL_LAPACK_FUNC( clatps , CLATPS )
#define F77_clatrd AOCL_LAPACK_FUNC( clatrd , CLATRD )
#define F77_clatrs AOCL_LAPACK_FUNC( clatrs , CLATRS )
#define F77_clatrz AOCL_LAPACK_FUNC( clatrz , CLATRZ )
#define F77_clatsqr AOCL_LAPACK_FUNC( clatsqr , CLATSQR )
#define F77_clatzm AOCL_LAPACK_FUNC( clatzm , CLATZM )
#define F77_claunhr_col_getrfnp2 F77_FUNC( claunhr_col_getrfnp2 , CLAUNHR_COL_GETRFNP2 )
#define F77_claunhr_col_getrfnp F77_FUNC( claunhr_col_getrfnp , CLAUNHR_COL_GETRFNP )
#define F77_clauu2 AOCL_LAPACK_FUNC( clauu2 , CLAUU2 )
#define F77_clauum AOCL_LAPACK_FUNC( clauum , CLAUUM )
#define F77_cla_wwaddw F77_FUNC( cla_wwaddw , CLA_WWADDW )
#define F77_cpbcon AOCL_LAPACK_FUNC( cpbcon , CPBCON )
#define F77_cpbequ AOCL_LAPACK_FUNC( cpbequ , CPBEQU )
#define F77_cpbrfs AOCL_LAPACK_FUNC( cpbrfs , CPBRFS )
#define F77_cpbstf AOCL_LAPACK_FUNC( cpbstf , CPBSTF )
#define F77_cpbsv AOCL_LAPACK_FUNC( cpbsv , CPBSV )
#define F77_cpbsvx AOCL_LAPACK_FUNC( cpbsvx , CPBSVX )
#define F77_cpbtf2 AOCL_LAPACK_FUNC( cpbtf2 , CPBTF2 )
#define F77_cpbtrf AOCL_LAPACK_FUNC( cpbtrf , CPBTRF )
#define F77_cpbtrs AOCL_LAPACK_FUNC( cpbtrs , CPBTRS )
#define F77_cpftrf AOCL_LAPACK_FUNC( cpftrf , CPFTRF )
#define F77_cpftri AOCL_LAPACK_FUNC( cpftri , CPFTRI )
#define F77_cpftrs AOCL_LAPACK_FUNC( cpftrs , CPFTRS )
#define F77_cpocon AOCL_LAPACK_FUNC( cpocon , CPOCON )
#define F77_cpoequb AOCL_LAPACK_FUNC( cpoequb , CPOEQUB )
#define F77_cpoequ AOCL_LAPACK_FUNC( cpoequ , CPOEQU )
#define F77_cporfs AOCL_LAPACK_FUNC( cporfs , CPORFS )
#define F77_cporfsx F77_FUNC( cporfsx , CPORFSX )
#define F77_cposv AOCL_LAPACK_FUNC( cposv , CPOSV )
#define F77_cposvx AOCL_LAPACK_FUNC( cposvx , CPOSVX )
#define F77_cposvxx F77_FUNC( cposvxx , CPOSVXX )
#define F77_cpotf2 AOCL_LAPACK_FUNC( cpotf2 , CPOTF2 )
#define F77_cpotrf2 AOCL_LAPACK_FUNC( cpotrf2 , CPOTRF2 )
#define F77_cpotrf AOCL_LAPACK_FUNC( cpotrf , CPOTRF )
#define F77_cpotri AOCL_LAPACK_FUNC( cpotri , CPOTRI )
#define F77_cpotrs AOCL_LAPACK_FUNC( cpotrs , CPOTRS )
#define F77_cppcon AOCL_LAPACK_FUNC( cppcon , CPPCON )
#define F77_cppequ AOCL_LAPACK_FUNC( cppequ , CPPEQU )
#define F77_cpprfs AOCL_LAPACK_FUNC( cpprfs , CPPRFS )
#define F77_cppsv AOCL_LAPACK_FUNC( cppsv , CPPSV )
#define F77_cppsvx AOCL_LAPACK_FUNC( cppsvx , CPPSVX )
#define F77_cpptrf AOCL_LAPACK_FUNC( cpptrf , CPPTRF )
#define F77_cpptri AOCL_LAPACK_FUNC( cpptri , CPPTRI )
#define F77_cpptrs AOCL_LAPACK_FUNC( cpptrs , CPPTRS )
#define F77_cpstf2 AOCL_LAPACK_FUNC( cpstf2 , CPSTF2 )
#define F77_cpstrf AOCL_LAPACK_FUNC( cpstrf , CPSTRF )
#define F77_cptcon AOCL_LAPACK_FUNC( cptcon , CPTCON )
#define F77_cpteqr AOCL_LAPACK_FUNC( cpteqr , CPTEQR )
#define F77_cptrfs AOCL_LAPACK_FUNC( cptrfs , CPTRFS )
#define F77_cptsv AOCL_LAPACK_FUNC( cptsv , CPTSV )
#define F77_cptsvx AOCL_LAPACK_FUNC( cptsvx , CPTSVX )
#define F77_cpttrf AOCL_LAPACK_FUNC( cpttrf , CPTTRF )
#define F77_cpttrs AOCL_LAPACK_FUNC( cpttrs , CPTTRS )
#define F77_cptts2 AOCL_LAPACK_FUNC( cptts2 , CPTTS2 )
#define F77_crot AOCL_LAPACK_FUNC( crot , CROT )
#define F77_cspcon AOCL_LAPACK_FUNC( cspcon , CSPCON )
#define F77_cspmv AOCL_LAPACK_FUNC( cspmv , CSPMV )
#define F77_cspr AOCL_LAPACK_FUNC( cspr , CSPR )
#define F77_csprfs AOCL_LAPACK_FUNC( csprfs , CSPRFS )
#define F77_cspsv AOCL_LAPACK_FUNC( cspsv , CSPSV )
#define F77_cspsvx AOCL_LAPACK_FUNC( cspsvx , CSPSVX )
#define F77_csptrf AOCL_LAPACK_FUNC( csptrf , CSPTRF )
#define F77_csptri AOCL_LAPACK_FUNC( csptri , CSPTRI )
#define F77_csptrs AOCL_LAPACK_FUNC( csptrs , CSPTRS )
#define F77_csrscl AOCL_LAPACK_FUNC( csrscl , CSRSCL )
#define F77_cstedc AOCL_LAPACK_FUNC( cstedc , CSTEDC )
#define F77_cstegr AOCL_LAPACK_FUNC( cstegr , CSTEGR )
#define F77_cstein AOCL_LAPACK_FUNC( cstein , CSTEIN )
#define F77_cstemr AOCL_LAPACK_FUNC( cstemr , CSTEMR )
#define F77_csteqr AOCL_LAPACK_FUNC( csteqr , CSTEQR )
#define F77_csycon_3 AOCL_LAPACK_FUNC( csycon_3 , CSYCON_3 )
#define F77_csycon AOCL_LAPACK_FUNC( csycon , CSYCON )
#define F77_csycon_rook AOCL_LAPACK_FUNC( csycon_rook , CSYCON_ROOK )
#define F77_csyconv AOCL_LAPACK_FUNC( csyconv , CSYCONV )
#define F77_csyconvf AOCL_LAPACK_FUNC( csyconvf , CSYCONVF )
#define F77_csyconvf_rook AOCL_LAPACK_FUNC( csyconvf_rook , CSYCONVF_ROOK )
#define F77_csyequb AOCL_LAPACK_FUNC( csyequb , CSYEQUB )
#define F77_csymv AOCL_LAPACK_FUNC( csymv , CSYMV )
#define F77_csyr AOCL_LAPACK_FUNC( csyr , CSYR )
#define F77_csyrfs AOCL_LAPACK_FUNC( csyrfs , CSYRFS )
#define F77_csyrfsx F77_FUNC( csyrfsx , CSYRFSX )
#define F77_csysv_aa_2stage AOCL_LAPACK_FUNC( csysv_aa_2stage , CSYSV_AA_2STAGE )
#define F77_csysv_aa AOCL_LAPACK_FUNC( csysv_aa , CSYSV_AA )
#define F77_csysv AOCL_LAPACK_FUNC( csysv , CSYSV )
#define F77_csysv_rk AOCL_LAPACK_FUNC( csysv_rk , CSYSV_RK )
#define F77_csysv_rook AOCL_LAPACK_FUNC( csysv_rook , CSYSV_ROOK )
#define F77_csysvx AOCL_LAPACK_FUNC( csysvx , CSYSVX )
#define F77_csysvxx F77_FUNC( csysvxx , CSYSVXX )
#define F77_csyswapr AOCL_LAPACK_FUNC( csyswapr , CSYSWAPR )
#define F77_csytf2 AOCL_LAPACK_FUNC( csytf2 , CSYTF2 )
#define F77_csytf2_rk AOCL_LAPACK_FUNC( csytf2_rk , CSYTF2_RK )
#define F77_csytf2_rook AOCL_LAPACK_FUNC( csytf2_rook , CSYTF2_ROOK )
#define F77_csytrf_aa_2stage AOCL_LAPACK_FUNC( csytrf_aa_2stage , CSYTRF_AA_2STAGE )
#define F77_csytrf_aa AOCL_LAPACK_FUNC( csytrf_aa , CSYTRF_AA )
#define F77_csytrf AOCL_LAPACK_FUNC( csytrf , CSYTRF )
#define F77_csytrf_rk AOCL_LAPACK_FUNC( csytrf_rk , CSYTRF_RK )
#define F77_csytrf_rook AOCL_LAPACK_FUNC( csytrf_rook , CSYTRF_ROOK )
#define F77_csytri2 AOCL_LAPACK_FUNC( csytri2 , CSYTRI2 )
#define F77_csytri2x AOCL_LAPACK_FUNC( csytri2x , CSYTRI2X )
#define F77_csytri_3 AOCL_LAPACK_FUNC( csytri_3 , CSYTRI_3 )
#define F77_csytri_3x AOCL_LAPACK_FUNC( csytri_3x , CSYTRI_3X )
#define F77_csytri AOCL_LAPACK_FUNC( csytri , CSYTRI )
#define F77_csytri_rook AOCL_LAPACK_FUNC( csytri_rook , CSYTRI_ROOK )
#define F77_csytrs2 AOCL_LAPACK_FUNC( csytrs2 , CSYTRS2 )
#define F77_csytrs_3 AOCL_LAPACK_FUNC( csytrs_3 , CSYTRS_3 )
#define F77_csytrs_aa_2stage AOCL_LAPACK_FUNC( csytrs_aa_2stage , CSYTRS_AA_2STAGE )
#define F77_csytrs_aa AOCL_LAPACK_FUNC( csytrs_aa , CSYTRS_AA )
#define F77_csytrs AOCL_LAPACK_FUNC( csytrs , CSYTRS )
#define F77_csytrs_rook AOCL_LAPACK_FUNC( csytrs_rook , CSYTRS_ROOK )
#define F77_ctbcon AOCL_LAPACK_FUNC( ctbcon , CTBCON )
#define F77_ctbrfs AOCL_LAPACK_FUNC( ctbrfs , CTBRFS )
#define F77_ctbtrs AOCL_LAPACK_FUNC( ctbtrs , CTBTRS )
#define F77_ctfsm AOCL_LAPACK_FUNC( ctfsm , CTFSM )
#define F77_ctftri AOCL_LAPACK_FUNC( ctftri , CTFTRI )
#define F77_ctfttp AOCL_LAPACK_FUNC( ctfttp , CTFTTP )
#define F77_ctfttr AOCL_LAPACK_FUNC( ctfttr , CTFTTR )
#define F77_ctgevc AOCL_LAPACK_FUNC( ctgevc , CTGEVC )
#define F77_ctgex2 AOCL_LAPACK_FUNC( ctgex2 , CTGEX2 )
#define F77_ctgexc AOCL_LAPACK_FUNC( ctgexc , CTGEXC )
#define F77_ctgsen AOCL_LAPACK_FUNC( ctgsen , CTGSEN )
#define F77_ctgsja AOCL_LAPACK_FUNC( ctgsja , CTGSJA )
#define F77_ctgsna AOCL_LAPACK_FUNC( ctgsna , CTGSNA )
#define F77_ctgsy2 AOCL_LAPACK_FUNC( ctgsy2 , CTGSY2 )
#define F77_ctgsyl AOCL_LAPACK_FUNC( ctgsyl , CTGSYL )
#define F77_ctpcon AOCL_LAPACK_FUNC( ctpcon , CTPCON )
#define F77_ctplqt2 AOCL_LAPACK_FUNC( ctplqt2 , CTPLQT2 )
#define F77_ctplqt AOCL_LAPACK_FUNC( ctplqt , CTPLQT )
#define F77_ctpmlqt AOCL_LAPACK_FUNC( ctpmlqt , CTPMLQT )
#define F77_ctpmqrt AOCL_LAPACK_FUNC( ctpmqrt , CTPMQRT )
#define F77_ctpqrt2 AOCL_LAPACK_FUNC( ctpqrt2 , CTPQRT2 )
#define F77_ctpqrt AOCL_LAPACK_FUNC( ctpqrt , CTPQRT )
#define F77_ctprfb AOCL_LAPACK_FUNC( ctprfb , CTPRFB )
#define F77_ctprfs AOCL_LAPACK_FUNC( ctprfs , CTPRFS )
#define F77_ctptri AOCL_LAPACK_FUNC( ctptri , CTPTRI )
#define F77_ctptrs AOCL_LAPACK_FUNC( ctptrs , CTPTRS )
#define F77_ctpttf AOCL_LAPACK_FUNC( ctpttf , CTPTTF )
#define F77_ctpttr AOCL_LAPACK_FUNC( ctpttr , CTPTTR )
#define F77_ctrcon AOCL_LAPACK_FUNC( ctrcon , CTRCON )
#define F77_ctrevc3 AOCL_LAPACK_FUNC( ctrevc3 , CTREVC3 )
#define F77_ctrevc AOCL_LAPACK_FUNC( ctrevc , CTREVC )
#define F77_ctrexc AOCL_LAPACK_FUNC( ctrexc , CTREXC )
#define F77_ctrrfs AOCL_LAPACK_FUNC( ctrrfs , CTRRFS )
#define F77_ctrsen AOCL_LAPACK_FUNC( ctrsen , CTRSEN )
#define F77_ctrsna AOCL_LAPACK_FUNC( ctrsna , CTRSNA )
#define F77_ctrsyl AOCL_LAPACK_FUNC( ctrsyl , CTRSYL )
#define F77_ctrti2 AOCL_LAPACK_FUNC( ctrti2 , CTRTI2 )
#define F77_ctrtri AOCL_LAPACK_FUNC( ctrtri , CTRTRI )
#define F77_ctrtrs AOCL_LAPACK_FUNC( ctrtrs , CTRTRS )
#define F77_ctrttf AOCL_LAPACK_FUNC( ctrttf , CTRTTF )
#define F77_ctrttp AOCL_LAPACK_FUNC( ctrttp , CTRTTP )
#define F77_ctzrqf AOCL_LAPACK_FUNC( ctzrqf , CTZRQF )
#define F77_ctzrzf AOCL_LAPACK_FUNC( ctzrzf , CTZRZF )
#define F77_cunbdb1 AOCL_LAPACK_FUNC( cunbdb1 , CUNBDB1 )
#define F77_cunbdb2 AOCL_LAPACK_FUNC( cunbdb2 , CUNBDB2 )
#define F77_cunbdb3 AOCL_LAPACK_FUNC( cunbdb3 , CUNBDB3 )
#define F77_cunbdb4 AOCL_LAPACK_FUNC( cunbdb4 , CUNBDB4 )
#define F77_cunbdb5 AOCL_LAPACK_FUNC( cunbdb5 , CUNBDB5 )
#define F77_cunbdb6 AOCL_LAPACK_FUNC( cunbdb6 , CUNBDB6 )
#define F77_cunbdb AOCL_LAPACK_FUNC( cunbdb , CUNBDB )
#define F77_cuncsd2by1 AOCL_LAPACK_FUNC( cuncsd2by1 , CUNCSD2BY1 )
#define F77_cuncsd AOCL_LAPACK_FUNC( cuncsd , CUNCSD )
#define F77_cung2l AOCL_LAPACK_FUNC( cung2l , CUNG2L )
#define F77_cung2r AOCL_LAPACK_FUNC( cung2r , CUNG2R )
#define F77_cungbr AOCL_LAPACK_FUNC( cungbr , CUNGBR )
#define F77_cunghr AOCL_LAPACK_FUNC( cunghr , CUNGHR )
#define F77_cungl2 AOCL_LAPACK_FUNC( cungl2 , CUNGL2 )
#define F77_cunglq AOCL_LAPACK_FUNC( cunglq , CUNGLQ )
#define F77_cungql AOCL_LAPACK_FUNC( cungql , CUNGQL )
#define F77_cungqr AOCL_LAPACK_FUNC( cungqr , CUNGQR )
#define F77_cungr2 AOCL_LAPACK_FUNC( cungr2 , CUNGR2 )
#define F77_cungrq AOCL_LAPACK_FUNC( cungrq , CUNGRQ )
#define F77_cungtr AOCL_LAPACK_FUNC( cungtr , CUNGTR )
#define F77_cungtsqr AOCL_LAPACK_FUNC( cungtsqr , CUNGTSQR )
#define F77_cunhr_col F77_FUNC( cunhr_col , CUNHR_COL )
#define F77_cunm22 AOCL_LAPACK_FUNC( cunm22 , CUNM22 )
#define F77_cunm2l AOCL_LAPACK_FUNC( cunm2l , CUNM2L )
#define F77_cunm2r AOCL_LAPACK_FUNC( cunm2r , CUNM2R )
#define F77_cunmbr AOCL_LAPACK_FUNC( cunmbr , CUNMBR )
#define F77_cunmhr AOCL_LAPACK_FUNC( cunmhr , CUNMHR )
#define F77_cunml2 AOCL_LAPACK_FUNC( cunml2 , CUNML2 )
#define F77_cunmlq AOCL_LAPACK_FUNC( cunmlq , CUNMLQ )
#define F77_cunmql AOCL_LAPACK_FUNC( cunmql , CUNMQL )
#define F77_cunmqr AOCL_LAPACK_FUNC( cunmqr , CUNMQR )
#define F77_cunmr2 AOCL_LAPACK_FUNC( cunmr2 , CUNMR2 )
#define F77_cunmr3 AOCL_LAPACK_FUNC( cunmr3 , CUNMR3 )
#define F77_cunmrq AOCL_LAPACK_FUNC( cunmrq , CUNMRQ )
#define F77_cunmrz AOCL_LAPACK_FUNC( cunmrz , CUNMRZ )
#define F77_cunmtr AOCL_LAPACK_FUNC( cunmtr , CUNMTR )
#define F77_cupgtr AOCL_LAPACK_FUNC( cupgtr , CUPGTR )
#define F77_cupmtr AOCL_LAPACK_FUNC( cupmtr , CUPMTR )
#define F77_dbbcsd AOCL_LAPACK_FUNC( dbbcsd , DBBCSD )
#define F77_dbdsdc AOCL_LAPACK_FUNC( dbdsdc , DBDSDC )
#define F77_dbdsqr AOCL_LAPACK_FUNC( dbdsqr , DBDSQR )
#define F77_dbdsvdx AOCL_LAPACK_FUNC( dbdsvdx , DBDSVDX )
#define F77_dcombssq F77_FUNC( dcombssq , DCOMBSSQ )
#define F77_ddisna AOCL_LAPACK_FUNC( ddisna , DDISNA )
#define F77_dgbbrd AOCL_LAPACK_FUNC( dgbbrd , DGBBRD )
#define F77_dgbcon AOCL_LAPACK_FUNC( dgbcon , DGBCON )
#define F77_dgbequb AOCL_LAPACK_FUNC( dgbequb , DGBEQUB )
#define F77_dgbequ AOCL_LAPACK_FUNC( dgbequ , DGBEQU )
#define F77_dgbrfs AOCL_LAPACK_FUNC( dgbrfs , DGBRFS )
#define F77_dgbrfsx F77_FUNC( dgbrfsx , DGBRFSX )
#define F77_dgbsv AOCL_LAPACK_FUNC( dgbsv , DGBSV )
#define F77_dgbsvx AOCL_LAPACK_FUNC( dgbsvx , DGBSVX )
#define F77_dgbsvxx F77_FUNC( dgbsvxx , DGBSVXX )
#define F77_dgbtf2 AOCL_LAPACK_FUNC( dgbtf2 , DGBTF2 )
#define F77_dgbtrf AOCL_LAPACK_FUNC( dgbtrf , DGBTRF )
#define F77_dgbtrs AOCL_LAPACK_FUNC( dgbtrs , DGBTRS )
#define F77_dgebak AOCL_LAPACK_FUNC( dgebak , DGEBAK )
#define F77_dgebal AOCL_LAPACK_FUNC( dgebal , DGEBAL )
#define F77_dgebd2 AOCL_LAPACK_FUNC( dgebd2 , DGEBD2 )
#define F77_dgebrd AOCL_LAPACK_FUNC( dgebrd , DGEBRD )
#define F77_dgecon AOCL_LAPACK_FUNC( dgecon , DGECON )
#define F77_dgeequb AOCL_LAPACK_FUNC( dgeequb , DGEEQUB )
#define F77_dgeequ AOCL_LAPACK_FUNC( dgeequ , DGEEQU )
#define F77_dgees AOCL_LAPACK_FUNC( dgees , DGEES )
#define F77_dgeesx AOCL_LAPACK_FUNC( dgeesx , DGEESX )
#define F77_dgeev AOCL_LAPACK_FUNC( dgeev , DGEEV )
#define F77_dgeevx AOCL_LAPACK_FUNC( dgeevx , DGEEVX )
#define F77_dgegs AOCL_LAPACK_FUNC( dgegs , DGEGS )
#define F77_dgegv AOCL_LAPACK_FUNC( dgegv , DGEGV )
#define F77_dgehd2 AOCL_LAPACK_FUNC( dgehd2 , DGEHD2 )
#define F77_dgehrd AOCL_LAPACK_FUNC( dgehrd , DGEHRD )
#define F77_dgejsv AOCL_LAPACK_FUNC( dgejsv , DGEJSV )
#define F77_dgelq2 AOCL_LAPACK_FUNC( dgelq2 , DGELQ2 )
#define F77_dgelq AOCL_LAPACK_FUNC( dgelq , DGELQ )
#define F77_dgelqf AOCL_LAPACK_FUNC( dgelqf , DGELQF )
#define F77_dgelqt3 AOCL_LAPACK_FUNC( dgelqt3 , DGELQT3 )
#define F77_dgelqt AOCL_LAPACK_FUNC( dgelqt , DGELQT )
#define F77_dgelsd AOCL_LAPACK_FUNC( dgelsd , DGELSD )
#define F77_dgels AOCL_LAPACK_FUNC( dgels , DGELS )
#define F77_dgelss AOCL_LAPACK_FUNC( dgelss , DGELSS )
#define F77_dgelsx AOCL_LAPACK_FUNC( dgelsx , DGELSX )
#define F77_dgelsy AOCL_LAPACK_FUNC( dgelsy , DGELSY )
#define F77_dgemlq AOCL_LAPACK_FUNC( dgemlq , DGEMLQ )
#define F77_dgemlqt AOCL_LAPACK_FUNC( dgemlqt , DGEMLQT )
#define F77_dgemqr AOCL_LAPACK_FUNC( dgemqr , DGEMQR )
#define F77_dgemqrt AOCL_LAPACK_FUNC( dgemqrt , DGEMQRT )
#define F77_dgeql2 AOCL_LAPACK_FUNC( dgeql2 , DGEQL2 )
#define F77_dgeqlf AOCL_LAPACK_FUNC( dgeqlf , DGEQLF )
#define F77_dgeqp3 AOCL_LAPACK_FUNC( dgeqp3 , DGEQP3 )
#define F77_dgeqpf AOCL_LAPACK_FUNC( dgeqpf , DGEQPF )
#define F77_dgeqr2 AOCL_LAPACK_FUNC( dgeqr2 , DGEQR2 )
#define F77_dgeqr2p AOCL_LAPACK_FUNC( dgeqr2p , DGEQR2P )
#define F77_dgeqr AOCL_LAPACK_FUNC( dgeqr , DGEQR )
#define F77_dgeqrf AOCL_LAPACK_FUNC( dgeqrf , DGEQRF )
#define F77_dgeqrfp AOCL_LAPACK_FUNC( dgeqrfp , DGEQRFP )
#define F77_dgeqrt2 AOCL_LAPACK_FUNC( dgeqrt2 , DGEQRT2 )
#define F77_dgeqrt3 AOCL_LAPACK_FUNC( dgeqrt3 , DGEQRT3 )
#define F77_dgeqrt AOCL_LAPACK_FUNC( dgeqrt , DGEQRT )
#define F77_dgerfs AOCL_LAPACK_FUNC( dgerfs , DGERFS )
#define F77_dgerfsx F77_FUNC( dgerfsx , DGERFSX )
#define F77_dgerq2 AOCL_LAPACK_FUNC( dgerq2 , DGERQ2 )
#define F77_dgerqf AOCL_LAPACK_FUNC( dgerqf , DGERQF )
#define F77_dgesc2 AOCL_LAPACK_FUNC( dgesc2 , DGESC2 )
#define F77_dgesdd AOCL_LAPACK_FUNC( dgesdd , DGESDD )
#define F77_dgesvd AOCL_LAPACK_FUNC( dgesvd , DGESVD )
#define F77_dgesvdq AOCL_LAPACK_FUNC( dgesvdq , DGESVDQ )
#define F77_dgesvdx AOCL_LAPACK_FUNC( dgesvdx , DGESVDX )
#define F77_dgesv AOCL_LAPACK_FUNC( dgesv , DGESV )
#define F77_dgesvj AOCL_LAPACK_FUNC( dgesvj , DGESVJ )
#define F77_dgesvx AOCL_LAPACK_FUNC( dgesvx , DGESVX )
#define F77_dgesvxx F77_FUNC( dgesvxx , DGESVXX )
#define F77_dgetc2 AOCL_LAPACK_FUNC( dgetc2 , DGETC2 )
#define F77_dgetf2 AOCL_LAPACK_FUNC( dgetf2 , DGETF2 )
#define F77_dgetrf2 AOCL_LAPACK_FUNC( dgetrf2 , DGETRF2 )
#define F77_dgetrf AOCL_LAPACK_FUNC( dgetrf , DGETRF )
#define F77_dgetri AOCL_LAPACK_FUNC( dgetri , DGETRI )
#define F77_dgetrs AOCL_LAPACK_FUNC( dgetrs , DGETRS )
#define F77_dgetsls AOCL_LAPACK_FUNC( dgetsls , DGETSLS )
#define F77_dggbak AOCL_LAPACK_FUNC( dggbak , DGGBAK )
#define F77_dggbal AOCL_LAPACK_FUNC( dggbal , DGGBAL )
#define F77_dgges3 AOCL_LAPACK_FUNC( dgges3 , DGGES3 )
#define F77_dgges AOCL_LAPACK_FUNC( dgges , DGGES )
#define F77_dggesx AOCL_LAPACK_FUNC( dggesx , DGGESX )
#define F77_dggev3 AOCL_LAPACK_FUNC( dggev3 , DGGEV3 )
#define F77_dggev AOCL_LAPACK_FUNC( dggev , DGGEV )
#define F77_dggevx AOCL_LAPACK_FUNC( dggevx , DGGEVX )
#define F77_dggglm AOCL_LAPACK_FUNC( dggglm , DGGGLM )
#define F77_dgghd3 AOCL_LAPACK_FUNC( dgghd3 , DGGHD3 )
#define F77_dgghrd AOCL_LAPACK_FUNC( dgghrd , DGGHRD )
#define F77_dgglse AOCL_LAPACK_FUNC( dgglse , DGGLSE )
#define F77_dggqrf AOCL_LAPACK_FUNC( dggqrf , DGGQRF )
#define F77_dggrqf AOCL_LAPACK_FUNC( dggrqf , DGGRQF )
#define F77_dggsvd3 AOCL_LAPACK_FUNC( dggsvd3 , DGGSVD3 )
#define F77_dggsvd AOCL_LAPACK_FUNC( dggsvd , DGGSVD )
#define F77_dggsvp3 AOCL_LAPACK_FUNC( dggsvp3 , DGGSVP3 )
#define F77_dggsvp AOCL_LAPACK_FUNC( dggsvp , DGGSVP )
#define F77_dgsvj0 AOCL_LAPACK_FUNC( dgsvj0 , DGSVJ0 )
#define F77_dgsvj1 AOCL_LAPACK_FUNC( dgsvj1 , DGSVJ1 )
#define F77_dgtcon AOCL_LAPACK_FUNC( dgtcon , DGTCON )
#define F77_dgtrfs AOCL_LAPACK_FUNC( dgtrfs , DGTRFS )
#define F77_dgtsv AOCL_LAPACK_FUNC( dgtsv , DGTSV )
#define F77_dgtsvx AOCL_LAPACK_FUNC( dgtsvx , DGTSVX )
#define F77_dgttrf AOCL_LAPACK_FUNC( dgttrf , DGTTRF )
#define F77_dgttrs AOCL_LAPACK_FUNC( dgttrs , DGTTRS )
#define F77_dgtts2 AOCL_LAPACK_FUNC( dgtts2 , DGTTS2 )
#define F77_dhgeqz AOCL_LAPACK_FUNC( dhgeqz , DHGEQZ )
#define F77_dhsein AOCL_LAPACK_FUNC( dhsein , DHSEIN )
#define F77_dhseqr AOCL_LAPACK_FUNC( dhseqr , DHSEQR )
#define F77_disnan F77_FUNC( disnan , DISNAN )
#define F77_dlabad F77_FUNC( dlabad , DLABAD )
#define F77_dlabrd AOCL_LAPACK_FUNC( dlabrd , DLABRD )
#define F77_dlacn2 AOCL_LAPACK_FUNC( dlacn2 , DLACN2 )
#define F77_dlacon AOCL_LAPACK_FUNC( dlacon , DLACON )
#define F77_dlacpy AOCL_LAPACK_FUNC( dlacpy , DLACPY )
#define F77_dladiv F77_FUNC( dladiv , DLADIV )
#define F77_dlae2 F77_FUNC( dlae2 , DLAE2 )
#define F77_dlaebz AOCL_LAPACK_FUNC( dlaebz , DLAEBZ )
#define F77_dlaed0 AOCL_LAPACK_FUNC( dlaed0 , DLAED0 )
#define F77_dlaed1 AOCL_LAPACK_FUNC( dlaed1 , DLAED1 )
#define F77_dlaed2 AOCL_LAPACK_FUNC( dlaed2 , DLAED2 )
#define F77_dlaed3 AOCL_LAPACK_FUNC( dlaed3 , DLAED3 )
#define F77_dlaed4 AOCL_LAPACK_FUNC( dlaed4 , DLAED4 )
#define F77_dlaed5 AOCL_LAPACK_FUNC( dlaed5 , DLAED5 )
#define F77_dlaed6 AOCL_LAPACK_FUNC( dlaed6 , DLAED6 )
#define F77_dlaed7 AOCL_LAPACK_FUNC( dlaed7 , DLAED7 )
#define F77_dlaed8 AOCL_LAPACK_FUNC( dlaed8 , DLAED8 )
#define F77_dlaed9 AOCL_LAPACK_FUNC( dlaed9 , DLAED9 )
#define F77_dlaeda AOCL_LAPACK_FUNC( dlaeda , DLAEDA )
#define F77_dlaein AOCL_LAPACK_FUNC( dlaein , DLAEIN )
#define F77_dlaev2 F77_FUNC( dlaev2 , DLAEV2 )
#define F77_dlaexc AOCL_LAPACK_FUNC( dlaexc , DLAEXC )
#define F77_dlag2 AOCL_LAPACK_FUNC( dlag2 , DLAG2 )
#define F77_dlag2s AOCL_LAPACK_FUNC( dlag2s , DLAG2S )
#define F77_dla_gbamv F77_FUNC( dla_gbamv , DLA_GBAMV )
#define F77_dla_gbrcond F77_FUNC( dla_gbrcond , DLA_GBRCOND )
#define F77_dla_gbrfsx_extended F77_FUNC( dla_gbrfsx_extended , DLA_GBRFSX_EXTENDED )
#define F77_dla_gbrpvgrw F77_FUNC( dla_gbrpvgrw , DLA_GBRPVGRW )
#define F77_dla_geamv F77_FUNC( dla_geamv , DLA_GEAMV )
#define F77_dla_gercond F77_FUNC( dla_gercond , DLA_GERCOND )
#define F77_dla_gerfsx_extended F77_FUNC( dla_gerfsx_extended , DLA_GERFSX_EXTENDED )
#define F77_dla_gerpvgrw F77_FUNC( dla_gerpvgrw , DLA_GERPVGRW )
#define F77_dlags2 F77_FUNC( dlags2 , DLAGS2 )
#define F77_dlagtf AOCL_LAPACK_FUNC( dlagtf , DLAGTF )
#define F77_dlagtm AOCL_LAPACK_FUNC( dlagtm , DLAGTM )
#define F77_dlagts AOCL_LAPACK_FUNC( dlagts , DLAGTS )
#define F77_dlagv2 AOCL_LAPACK_FUNC( dlagv2 , DLAGV2 )
#define F77_dlahqr AOCL_LAPACK_FUNC( dlahqr , DLAHQR )
#define F77_dlahr2 AOCL_LAPACK_FUNC( dlahr2 , DLAHR2 )
#define F77_dlahrd AOCL_LAPACK_FUNC( dlahrd , DLAHRD )
#define F77_dlaic1 AOCL_LAPACK_FUNC( dlaic1 , DLAIC1 )
#define F77_dlaisnan F77_FUNC( dlaisnan , DLAISNAN )
#define F77_dla_lin_berr F77_FUNC( dla_lin_berr , DLA_LIN_BERR )
#define F77_dlaln2 AOCL_LAPACK_FUNC( dlaln2 , DLALN2 )
#define F77_dlals0 AOCL_LAPACK_FUNC( dlals0 , DLALS0 )
#define F77_dlalsa AOCL_LAPACK_FUNC( dlalsa , DLALSA )
#define F77_dlalsd AOCL_LAPACK_FUNC( dlalsd , DLALSD )
#define F77_dlamrg AOCL_LAPACK_FUNC( dlamrg , DLAMRG )
#define F77_dlamswlq AOCL_LAPACK_FUNC( dlamswlq , DLAMSWLQ )
#define F77_dlamtsqr AOCL_LAPACK_FUNC( dlamtsqr , DLAMTSQR )
#define F77_dlaneg AOCL_LAPACK_FUNC( dlaneg , DLANEG )
#define F77_dlangb AOCL_LAPACK_FUNC( dlangb , DLANGB )
#define F77_dlange AOCL_LAPACK_FUNC( dlange , DLANGE )
#define F77_dlangt AOCL_LAPACK_FUNC( dlangt , DLANGT )
#define F77_dlanhs AOCL_LAPACK_FUNC( dlanhs , DLANHS )
#define F77_dlansb AOCL_LAPACK_FUNC( dlansb , DLANSB )
#define F77_dlansf AOCL_LAPACK_FUNC( dlansf , DLANSF )
#define F77_dlansp AOCL_LAPACK_FUNC( dlansp , DLANSP )
#define F77_dlanst AOCL_LAPACK_FUNC( dlanst , DLANST )
#define F77_dlansy AOCL_LAPACK_FUNC( dlansy , DLANSY )
#define F77_dlantb AOCL_LAPACK_FUNC( dlantb , DLANTB )
#define F77_dlantp AOCL_LAPACK_FUNC( dlantp , DLANTP )
#define F77_dlantr AOCL_LAPACK_FUNC( dlantr , DLANTR )
#define F77_dlanv2 F77_FUNC( dlanv2 , DLANV2 )
#define F77_dlaorhr_col_getrfnp2 F77_FUNC( dlaorhr_col_getrfnp2 , DLAORHR_COL_GETRFNP2 )
#define F77_dlaorhr_col_getrfnp F77_FUNC( dlaorhr_col_getrfnp , DLAORHR_COL_GETRFNP )
#define F77_dlapll AOCL_LAPACK_FUNC( dlapll , DLAPLL )
#define F77_dlapmr AOCL_LAPACK_FUNC( dlapmr , DLAPMR )
#define F77_dlapmt AOCL_LAPACK_FUNC( dlapmt , DLAPMT )
#define F77_dla_porcond F77_FUNC( dla_porcond , DLA_PORCOND )
#define F77_dla_porfsx_extended F77_FUNC( dla_porfsx_extended , DLA_PORFSX_EXTENDED )
#define F77_dla_porpvgrw F77_FUNC( dla_porpvgrw , DLA_PORPVGRW )
#define F77_dlapy2 F77_FUNC( dlapy2 , DLAPY2 )
#define F77_dlapy3 F77_FUNC( dlapy3 , DLAPY3 )
#define F77_dlaqgb AOCL_LAPACK_FUNC( dlaqgb , DLAQGB )
#define F77_dlaqge AOCL_LAPACK_FUNC( dlaqge , DLAQGE )
#define F77_dlaqp2 AOCL_LAPACK_FUNC( dlaqp2 , DLAQP2 )
#define F77_dlaqps AOCL_LAPACK_FUNC( dlaqps , DLAQPS )
#define F77_dlaqr0 AOCL_LAPACK_FUNC( dlaqr0 , DLAQR0 )
#define F77_dlaqr1 AOCL_LAPACK_FUNC( dlaqr1 , DLAQR1 )
#define F77_dlaqr2 AOCL_LAPACK_FUNC( dlaqr2 , DLAQR2 )
#define F77_dlaqr3 AOCL_LAPACK_FUNC( dlaqr3 , DLAQR3 )
#define F77_dlaqr4 AOCL_LAPACK_FUNC( dlaqr4 , DLAQR4 )
#define F77_dlaqr5 AOCL_LAPACK_FUNC( dlaqr5 , DLAQR5 )
#define F77_dlaqsb AOCL_LAPACK_FUNC( dlaqsb , DLAQSB )
#define F77_dlaqsp AOCL_LAPACK_FUNC( dlaqsp , DLAQSP )
#define F77_dlaqsy AOCL_LAPACK_FUNC( dlaqsy , DLAQSY )
#define F77_dlaqtr AOCL_LAPACK_FUNC( dlaqtr , DLAQTR )
#define F77_dlar1v AOCL_LAPACK_FUNC( dlar1v , DLAR1V )
#define F77_dlar2v AOCL_LAPACK_FUNC( dlar2v , DLAR2V )
#define F77_dlarfb AOCL_LAPACK_FUNC( dlarfb , DLARFB )
#define F77_dlarf AOCL_LAPACK_FUNC( dlarf , DLARF )
#define F77_dlarfg AOCL_LAPACK_FUNC( dlarfg , DLARFG )
#define F77_dlarfgp AOCL_LAPACK_FUNC( dlarfgp , DLARFGP )
#define F77_dlarft AOCL_LAPACK_FUNC( dlarft , DLARFT )
#define F77_dlarfx AOCL_LAPACK_FUNC( dlarfx , DLARFX )
#define F77_dlarfy AOCL_LAPACK_FUNC( dlarfy , DLARFY )
#define F77_dlargv AOCL_LAPACK_FUNC( dlargv , DLARGV )
#define F77_dlarnv AOCL_LAPACK_FUNC( dlarnv , DLARNV )
#define F77_dlarra AOCL_LAPACK_FUNC( dlarra , DLARRA )
#define F77_dlarrb AOCL_LAPACK_FUNC( dlarrb , DLARRB )
#define F77_dlarrc AOCL_LAPACK_FUNC( dlarrc , DLARRC )
#define F77_dlarrd AOCL_LAPACK_FUNC( dlarrd , DLARRD )
#define F77_dlarre AOCL_LAPACK_FUNC( dlarre , DLARRE )
#define F77_dlarrf AOCL_LAPACK_FUNC( dlarrf , DLARRF )
#define F77_dlarrj AOCL_LAPACK_FUNC( dlarrj , DLARRJ )
#define F77_dlarrk AOCL_LAPACK_FUNC( dlarrk , DLARRK )
#define F77_dlarrr AOCL_LAPACK_FUNC( dlarrr , DLARRR )
#define F77_dlarrv AOCL_LAPACK_FUNC( dlarrv , DLARRV )
#define F77_dlarscl2 AOCL_LAPACK_FUNC( dlarscl2 , DLARSCL2 )
#define F77_dlartg F77_FUNC( dlartg , DLARTG )
#define F77_dlartgp F77_FUNC( dlartgp , DLARTGP )
#define F77_dlartgs F77_FUNC( dlartgs , DLARTGS )
#define F77_dlartv AOCL_LAPACK_FUNC( dlartv , DLARTV )
#define F77_dlaruv AOCL_LAPACK_FUNC( dlaruv , DLARUV )
#define F77_dlarzb AOCL_LAPACK_FUNC( dlarzb , DLARZB )
#define F77_dlarz AOCL_LAPACK_FUNC( dlarz , DLARZ )
#define F77_dlarzt AOCL_LAPACK_FUNC( dlarzt , DLARZT )
#define F77_dlas2 F77_FUNC( dlas2 , DLAS2 )
#define F77_dlascl2 AOCL_LAPACK_FUNC( dlascl2 , DLASCL2 )
#define F77_dlascl AOCL_LAPACK_FUNC( dlascl , DLASCL )
#define F77_dlasd0 AOCL_LAPACK_FUNC( dlasd0 , DLASD0 )
#define F77_dlasd1 AOCL_LAPACK_FUNC( dlasd1 , DLASD1 )
#define F77_dlasd2 AOCL_LAPACK_FUNC( dlasd2 , DLASD2 )
#define F77_dlasd3 AOCL_LAPACK_FUNC( dlasd3 , DLASD3 )
#define F77_dlasd4 AOCL_LAPACK_FUNC( dlasd4 , DLASD4 )
#define F77_dlasd5 AOCL_LAPACK_FUNC( dlasd5 , DLASD5 )
#define F77_dlasd6 AOCL_LAPACK_FUNC( dlasd6 , DLASD6 )
#define F77_dlasd7 AOCL_LAPACK_FUNC( dlasd7 , DLASD7 )
#define F77_dlasd8 AOCL_LAPACK_FUNC( dlasd8 , DLASD8 )
#define F77_dlasda AOCL_LAPACK_FUNC( dlasda , DLASDA )
#define F77_dlasdq AOCL_LAPACK_FUNC( dlasdq , DLASDQ )
#define F77_dlasdt AOCL_LAPACK_FUNC( dlasdt , DLASDT )
#define F77_dlaset AOCL_LAPACK_FUNC( dlaset , DLASET )
#define F77_dlasq1 AOCL_LAPACK_FUNC( dlasq1 , DLASQ1 )
#define F77_dlasq2 AOCL_LAPACK_FUNC( dlasq2 , DLASQ2 )
#define F77_dlasq3 AOCL_LAPACK_FUNC( dlasq3 , DLASQ3 )
#define F77_dlasq4 AOCL_LAPACK_FUNC( dlasq4 , DLASQ4 )
#define F77_dlasq5 AOCL_LAPACK_FUNC( dlasq5 , DLASQ5 )
#define F77_dlasq6 AOCL_LAPACK_FUNC( dlasq6 , DLASQ6 )
#define F77_dlasr AOCL_LAPACK_FUNC( dlasr , DLASR )
#define F77_dlasrt AOCL_LAPACK_FUNC( dlasrt , DLASRT )
#define F77_dlassq AOCL_LAPACK_FUNC( dlassq , DLASSQ )
#define F77_dlasv2 F77_FUNC( dlasv2 , DLASV2 )
#define F77_dlaswlq AOCL_LAPACK_FUNC( dlaswlq , DLASWLQ )
#define F77_dlaswp AOCL_LAPACK_FUNC( dlaswp , DLASWP )
#define F77_dlasy2 AOCL_LAPACK_FUNC( dlasy2 , DLASY2 )
#define F77_dla_syamv F77_FUNC( dla_syamv , DLA_SYAMV )
#define F77_dlasyf_aa AOCL_LAPACK_FUNC( dlasyf_aa , DLASYF_AA )
#define F77_dlasyf AOCL_LAPACK_FUNC( dlasyf , DLASYF )
#define F77_dlasyf_rk AOCL_LAPACK_FUNC( dlasyf_rk , DLASYF_RK )
#define F77_dlasyf_rook AOCL_LAPACK_FUNC( dlasyf_rook , DLASYF_ROOK )
#define F77_dla_syrcond F77_FUNC( dla_syrcond , DLA_SYRCOND )
#define F77_dla_syrfsx_extended F77_FUNC( dla_syrfsx_extended , DLA_SYRFSX_EXTENDED )
#define F77_dla_syrpvgrw F77_FUNC( dla_syrpvgrw , DLA_SYRPVGRW )
#define F77_dlat2s AOCL_LAPACK_FUNC( dlat2s , DLAT2S )
#define F77_dlatbs AOCL_LAPACK_FUNC( dlatbs , DLATBS )
#define F77_dlatdf AOCL_LAPACK_FUNC( dlatdf , DLATDF )
#define F77_dlatps AOCL_LAPACK_FUNC( dlatps , DLATPS )
#define F77_dlatrd AOCL_LAPACK_FUNC( dlatrd , DLATRD )
#define F77_dlatrs AOCL_LAPACK_FUNC( dlatrs , DLATRS )
#define F77_dlatrz AOCL_LAPACK_FUNC( dlatrz , DLATRZ )
#define F77_dlatsqr AOCL_LAPACK_FUNC( dlatsqr , DLATSQR )
#define F77_dlatzm AOCL_LAPACK_FUNC( dlatzm , DLATZM )
#define F77_dlauu2 AOCL_LAPACK_FUNC( dlauu2 , DLAUU2 )
#define F77_dlauum AOCL_LAPACK_FUNC( dlauum , DLAUUM )
#define F77_dla_wwaddw F77_FUNC( dla_wwaddw , DLA_WWADDW )
#define F77_dopmtr AOCL_LAPACK_FUNC( dopmtr , DOPMTR )
#define F77_dorbdb1 AOCL_LAPACK_FUNC( dorbdb1 , DORBDB1 )
#define F77_dorbdb2 AOCL_LAPACK_FUNC( dorbdb2 , DORBDB2 )
#define F77_dorbdb3 AOCL_LAPACK_FUNC( dorbdb3 , DORBDB3 )
#define F77_dorbdb4 AOCL_LAPACK_FUNC( dorbdb4 , DORBDB4 )
#define F77_dorbdb5 AOCL_LAPACK_FUNC( dorbdb5 , DORBDB5 )
#define F77_dorbdb6 AOCL_LAPACK_FUNC( dorbdb6 , DORBDB6 )
#define F77_dorbdb AOCL_LAPACK_FUNC( dorbdb , DORBDB )
#define F77_dorg2l AOCL_LAPACK_FUNC( dorg2l , DORG2L )
#define F77_dorg2r AOCL_LAPACK_FUNC( dorg2r , DORG2R )
#define F77_dorgbr AOCL_LAPACK_FUNC( dorgbr , DORGBR )
#define F77_dorgl2 AOCL_LAPACK_FUNC( dorgl2 , DORGL2 )
#define F77_dorglq AOCL_LAPACK_FUNC( dorglq , DORGLQ )
#define F77_dorgql AOCL_LAPACK_FUNC( dorgql , DORGQL )
#define F77_dorgqr AOCL_LAPACK_FUNC( dorgqr , DORGQR )
#define F77_dorgr2 AOCL_LAPACK_FUNC( dorgr2 , DORGR2 )
#define F77_dorgrq AOCL_LAPACK_FUNC( dorgrq , DORGRQ )
#define F77_dorgtr AOCL_LAPACK_FUNC( dorgtr , DORGTR )
#define F77_dorgtsqr AOCL_LAPACK_FUNC( dorgtsqr , DORGTSQR )
#define F77_dorhr_col F77_FUNC( dorhr_col , DORHR_COL )
#define F77_dorm22 AOCL_LAPACK_FUNC( dorm22 , DORM22 )
#define F77_dorm2l AOCL_LAPACK_FUNC( dorm2l , DORM2L )
#define F77_dorm2r AOCL_LAPACK_FUNC( dorm2r , DORM2R )
#define F77_dormbr AOCL_LAPACK_FUNC( dormbr , DORMBR )
#define F77_dorml2 AOCL_LAPACK_FUNC( dorml2 , DORML2 )
#define F77_dormlq AOCL_LAPACK_FUNC( dormlq , DORMLQ )
#define F77_dormql AOCL_LAPACK_FUNC( dormql , DORMQL )
#define F77_dormqr AOCL_LAPACK_FUNC( dormqr , DORMQR )
#define F77_dormr2 AOCL_LAPACK_FUNC( dormr2 , DORMR2 )
#define F77_dormr3 AOCL_LAPACK_FUNC( dormr3 , DORMR3 )
#define F77_dormrq AOCL_LAPACK_FUNC( dormrq , DORMRQ )
#define F77_dormrz AOCL_LAPACK_FUNC( dormrz , DORMRZ )
#define F77_dormtr AOCL_LAPACK_FUNC( dormtr , DORMTR )
#define F77_dpbcon AOCL_LAPACK_FUNC( dpbcon , DPBCON )
#define F77_dpbequ AOCL_LAPACK_FUNC( dpbequ , DPBEQU )
#define F77_dpbrfs AOCL_LAPACK_FUNC( dpbrfs , DPBRFS )
#define F77_dpbstf AOCL_LAPACK_FUNC( dpbstf , DPBSTF )
#define F77_dpbsv AOCL_LAPACK_FUNC( dpbsv , DPBSV )
#define F77_dpbsvx AOCL_LAPACK_FUNC( dpbsvx , DPBSVX )
#define F77_dpbtf2 AOCL_LAPACK_FUNC( dpbtf2 , DPBTF2 )
#define F77_dpbtrf AOCL_LAPACK_FUNC( dpbtrf , DPBTRF )
#define F77_dpbtrs AOCL_LAPACK_FUNC( dpbtrs , DPBTRS )
#define F77_dpftrf AOCL_LAPACK_FUNC( dpftrf , DPFTRF )
#define F77_dpftri AOCL_LAPACK_FUNC( dpftri , DPFTRI )
#define F77_dpftrs AOCL_LAPACK_FUNC( dpftrs , DPFTRS )
#define F77_dpocon AOCL_LAPACK_FUNC( dpocon , DPOCON )
#define F77_dpoequb AOCL_LAPACK_FUNC( dpoequb , DPOEQUB )
#define F77_dpoequ AOCL_LAPACK_FUNC( dpoequ , DPOEQU )
#define F77_dporfs AOCL_LAPACK_FUNC( dporfs , DPORFS )
#define F77_dporfsx F77_FUNC( dporfsx , DPORFSX )
#define F77_dposv AOCL_LAPACK_FUNC( dposv , DPOSV )
#define F77_dposvx AOCL_LAPACK_FUNC( dposvx , DPOSVX )
#define F77_dposvxx F77_FUNC( dposvxx , DPOSVXX )
#define F77_dpotf2 AOCL_LAPACK_FUNC( dpotf2 , DPOTF2 )
#define F77_dpotrf2 AOCL_LAPACK_FUNC( dpotrf2 , DPOTRF2 )
#define F77_dpotrf AOCL_LAPACK_FUNC( dpotrf , DPOTRF )
#define F77_dpotri AOCL_LAPACK_FUNC( dpotri , DPOTRI )
#define F77_dpotrs AOCL_LAPACK_FUNC( dpotrs , DPOTRS )
#define F77_dppcon AOCL_LAPACK_FUNC( dppcon , DPPCON )
#define F77_dppequ AOCL_LAPACK_FUNC( dppequ , DPPEQU )
#define F77_dpprfs AOCL_LAPACK_FUNC( dpprfs , DPPRFS )
#define F77_dppsv AOCL_LAPACK_FUNC( dppsv , DPPSV )
#define F77_dppsvx AOCL_LAPACK_FUNC( dppsvx , DPPSVX )
#define F77_dpptrf AOCL_LAPACK_FUNC( dpptrf , DPPTRF )
#define F77_dpptri AOCL_LAPACK_FUNC( dpptri , DPPTRI )
#define F77_dpptrs AOCL_LAPACK_FUNC( dpptrs , DPPTRS )
#define F77_dpstf2 AOCL_LAPACK_FUNC( dpstf2 , DPSTF2 )
#define F77_dpstrf AOCL_LAPACK_FUNC( dpstrf , DPSTRF )
#define F77_dptcon AOCL_LAPACK_FUNC( dptcon , DPTCON )
#define F77_dpteqr AOCL_LAPACK_FUNC( dpteqr , DPTEQR )
#define F77_dptrfs AOCL_LAPACK_FUNC( dptrfs , DPTRFS )
#define F77_dptsv AOCL_LAPACK_FUNC( dptsv , DPTSV )
#define F77_dptsvx AOCL_LAPACK_FUNC( dptsvx , DPTSVX )
#define F77_dpttrf AOCL_LAPACK_FUNC( dpttrf , DPTTRF )
#define F77_dpttrs AOCL_LAPACK_FUNC( dpttrs , DPTTRS )
#define F77_dptts2 AOCL_LAPACK_FUNC( dptts2 , DPTTS2 )
#define F77_drscl AOCL_LAPACK_FUNC( drscl , DRSCL )
#define F77_dsb2st_kernels F77_FUNC( dsb2st_kernels , DSB2ST_KERNELS )
#define F77_dsbev_2stage AOCL_LAPACK_FUNC( dsbev_2stage , DSBEV_2STAGE )
#define F77_dsbevd_2stage AOCL_LAPACK_FUNC( dsbevd_2stage , DSBEVD_2STAGE )
#define F77_dsbevd AOCL_LAPACK_FUNC( dsbevd , DSBEVD )
#define F77_dsbev AOCL_LAPACK_FUNC( dsbev , DSBEV )
#define F77_dsbevx_2stage AOCL_LAPACK_FUNC( dsbevx_2stage , DSBEVX_2STAGE )
#define F77_dsbevx AOCL_LAPACK_FUNC( dsbevx , DSBEVX )
#define F77_dsbgst AOCL_LAPACK_FUNC( dsbgst , DSBGST )
#define F77_dsbgvd AOCL_LAPACK_FUNC( dsbgvd , DSBGVD )
#define F77_dsbgv AOCL_LAPACK_FUNC( dsbgv , DSBGV )
#define F77_dsbgvx AOCL_LAPACK_FUNC( dsbgvx , DSBGVX )
#define F77_dsbtrd AOCL_LAPACK_FUNC( dsbtrd , DSBTRD )
#define F77_dsfrk AOCL_LAPACK_FUNC( dsfrk , DSFRK )
#define F77_dsgesv AOCL_LAPACK_FUNC( dsgesv , DSGESV )
#define F77_dspcon AOCL_LAPACK_FUNC( dspcon , DSPCON )
#define F77_dspevd AOCL_LAPACK_FUNC( dspevd , DSPEVD )
#define F77_dspev AOCL_LAPACK_FUNC( dspev , DSPEV )
#define F77_dspevx AOCL_LAPACK_FUNC( dspevx , DSPEVX )
#define F77_dspgst AOCL_LAPACK_FUNC( dspgst , DSPGST )
#define F77_dspgvd AOCL_LAPACK_FUNC( dspgvd , DSPGVD )
#define F77_dspgv AOCL_LAPACK_FUNC( dspgv , DSPGV )
#define F77_dspgvx AOCL_LAPACK_FUNC( dspgvx , DSPGVX )
#define F77_dsposv AOCL_LAPACK_FUNC( dsposv , DSPOSV )
#define F77_dsprfs AOCL_LAPACK_FUNC( dsprfs , DSPRFS )
#define F77_dspsv AOCL_LAPACK_FUNC( dspsv , DSPSV )
#define F77_dspsvx AOCL_LAPACK_FUNC( dspsvx , DSPSVX )
#define F77_dsptrd AOCL_LAPACK_FUNC( dsptrd , DSPTRD )
#define F77_dsptrf AOCL_LAPACK_FUNC( dsptrf , DSPTRF )
#define F77_dsptri AOCL_LAPACK_FUNC( dsptri , DSPTRI )
#define F77_dsptrs AOCL_LAPACK_FUNC( dsptrs , DSPTRS )
#define F77_dstebz AOCL_LAPACK_FUNC( dstebz , DSTEBZ )
#define F77_dstedc AOCL_LAPACK_FUNC( dstedc , DSTEDC )
#define F77_dstegr AOCL_LAPACK_FUNC( dstegr , DSTEGR )
#define F77_dstein AOCL_LAPACK_FUNC( dstein , DSTEIN )
#define F77_dstemr AOCL_LAPACK_FUNC( dstemr , DSTEMR )
#define F77_dsteqr AOCL_LAPACK_FUNC( dsteqr , DSTEQR )
#define F77_dsterf AOCL_LAPACK_FUNC( dsterf , DSTERF )
#define F77_dstevd AOCL_LAPACK_FUNC( dstevd , DSTEVD )
#define F77_dstev AOCL_LAPACK_FUNC( dstev , DSTEV )
#define F77_dstevr AOCL_LAPACK_FUNC( dstevr , DSTEVR )
#define F77_dstevx AOCL_LAPACK_FUNC( dstevx , DSTEVX )
#define F77_dsycon_3 AOCL_LAPACK_FUNC( dsycon_3 , DSYCON_3 )
#define F77_dsycon AOCL_LAPACK_FUNC( dsycon , DSYCON )
#define F77_dsycon_rook AOCL_LAPACK_FUNC( dsycon_rook , DSYCON_ROOK )
#define F77_dsyconv AOCL_LAPACK_FUNC( dsyconv , DSYCONV )
#define F77_dsyconvf AOCL_LAPACK_FUNC( dsyconvf , DSYCONVF )
#define F77_dsyconvf_rook AOCL_LAPACK_FUNC( dsyconvf_rook , DSYCONVF_ROOK )
#define F77_dsyequb AOCL_LAPACK_FUNC( dsyequb , DSYEQUB )
#define F77_dsyev_2stage AOCL_LAPACK_FUNC( dsyev_2stage , DSYEV_2STAGE )
#define F77_dsyevd_2stage AOCL_LAPACK_FUNC( dsyevd_2stage , DSYEVD_2STAGE )
#define F77_dsyevd AOCL_LAPACK_FUNC( dsyevd , DSYEVD )
#define F77_dsyev  AOCL_LAPACK_FUNC( dsyev  , DSYEV  )
#define F77_dsyev AOCL_LAPACK_FUNC( dsyev , DSYEV )
#define F77_dsyevr_2stage AOCL_LAPACK_FUNC( dsyevr_2stage , DSYEVR_2STAGE )
#define F77_dsyevr AOCL_LAPACK_FUNC( dsyevr , DSYEVR )
#define F77_dsyevx_2stage AOCL_LAPACK_FUNC( dsyevx_2stage , DSYEVX_2STAGE )
#define F77_dsyevx AOCL_LAPACK_FUNC( dsyevx , DSYEVX )
#define F77_dsygs2 AOCL_LAPACK_FUNC( dsygs2 , DSYGS2 )
#define F77_dsygst AOCL_LAPACK_FUNC( dsygst , DSYGST )
#define F77_dsygv_2stage AOCL_LAPACK_FUNC( dsygv_2stage , DSYGV_2STAGE )
#define F77_dsygvd AOCL_LAPACK_FUNC( dsygvd , DSYGVD )
#define F77_dsygv AOCL_LAPACK_FUNC( dsygv , DSYGV )
#define F77_dsygvx AOCL_LAPACK_FUNC( dsygvx , DSYGVX )
#define F77_dsyrfs AOCL_LAPACK_FUNC( dsyrfs , DSYRFS )
#define F77_dsyrfsx F77_FUNC( dsyrfsx , DSYRFSX )
#define F77_dsysv_aa_2stage AOCL_LAPACK_FUNC( dsysv_aa_2stage , DSYSV_AA_2STAGE )
#define F77_dsysv_aa AOCL_LAPACK_FUNC( dsysv_aa , DSYSV_AA )
#define F77_dsysv AOCL_LAPACK_FUNC( dsysv , DSYSV )
#define F77_dsysv_rk AOCL_LAPACK_FUNC( dsysv_rk , DSYSV_RK )
#define F77_dsysv_rook AOCL_LAPACK_FUNC( dsysv_rook , DSYSV_ROOK )
#define F77_dsysvx AOCL_LAPACK_FUNC( dsysvx , DSYSVX )
#define F77_dsysvxx F77_FUNC( dsysvxx , DSYSVXX )
#define F77_dsyswapr AOCL_LAPACK_FUNC( dsyswapr , DSYSWAPR )
#define F77_dsytd2 AOCL_LAPACK_FUNC( dsytd2 , DSYTD2 )
#define F77_dsytf2 AOCL_LAPACK_FUNC( dsytf2 , DSYTF2 )
#define F77_dsytf2_rk AOCL_LAPACK_FUNC( dsytf2_rk , DSYTF2_RK )
#define F77_dsytf2_rook AOCL_LAPACK_FUNC( dsytf2_rook , DSYTF2_ROOK )
#define F77_dsytrd_2stage AOCL_LAPACK_FUNC( dsytrd_2stage , DSYTRD_2STAGE )
#define F77_dsytrd AOCL_LAPACK_FUNC( dsytrd , DSYTRD )
#define F77_dsytrd_sb2st AOCL_LAPACK_FUNC( dsytrd_sb2st , DSYTRD_SB2ST )
#define F77_dsytrd_sy2sb AOCL_LAPACK_FUNC( dsytrd_sy2sb , DSYTRD_SY2SB )
#define F77_dsytrf_aa_2stage AOCL_LAPACK_FUNC( dsytrf_aa_2stage , DSYTRF_AA_2STAGE )
#define F77_dsytrf_aa AOCL_LAPACK_FUNC( dsytrf_aa , DSYTRF_AA )
#define F77_dsytrf AOCL_LAPACK_FUNC( dsytrf , DSYTRF )
#define F77_dsytrf_rk AOCL_LAPACK_FUNC( dsytrf_rk , DSYTRF_RK )
#define F77_dsytrf_rook AOCL_LAPACK_FUNC( dsytrf_rook , DSYTRF_ROOK )
#define F77_dsytri2 AOCL_LAPACK_FUNC( dsytri2 , DSYTRI2 )
#define F77_dsytri2x AOCL_LAPACK_FUNC( dsytri2x , DSYTRI2X )
#define F77_dsytri_3 AOCL_LAPACK_FUNC( dsytri_3 , DSYTRI_3 )
#define F77_dsytri_3x AOCL_LAPACK_FUNC( dsytri_3x , DSYTRI_3X )
#define F77_dsytri AOCL_LAPACK_FUNC( dsytri , DSYTRI )
#define F77_dsytri_rook AOCL_LAPACK_FUNC( dsytri_rook , DSYTRI_ROOK )
#define F77_dsytrs2 AOCL_LAPACK_FUNC( dsytrs2 , DSYTRS2 )
#define F77_dsytrs_3 AOCL_LAPACK_FUNC( dsytrs_3 , DSYTRS_3 )
#define F77_dsytrs_aa_2stage AOCL_LAPACK_FUNC( dsytrs_aa_2stage , DSYTRS_AA_2STAGE )
#define F77_dsytrs_aa AOCL_LAPACK_FUNC( dsytrs_aa , DSYTRS_AA )
#define F77_dsytrs AOCL_LAPACK_FUNC( dsytrs , DSYTRS )
#define F77_dsytrs_rook AOCL_LAPACK_FUNC( dsytrs_rook , DSYTRS_ROOK )
#define F77_dtbcon AOCL_LAPACK_FUNC( dtbcon , DTBCON )
#define F77_dtbrfs AOCL_LAPACK_FUNC( dtbrfs , DTBRFS )
#define F77_dtbtrs AOCL_LAPACK_FUNC( dtbtrs , DTBTRS )
#define F77_dtfsm AOCL_LAPACK_FUNC( dtfsm , DTFSM )
#define F77_dtftri AOCL_LAPACK_FUNC( dtftri , DTFTRI )
#define F77_dtfttp AOCL_LAPACK_FUNC( dtfttp , DTFTTP )
#define F77_dtfttr AOCL_LAPACK_FUNC( dtfttr , DTFTTR )
#define F77_dtgevc AOCL_LAPACK_FUNC( dtgevc , DTGEVC )
#define F77_dtgex2 AOCL_LAPACK_FUNC( dtgex2 , DTGEX2 )
#define F77_dtgexc AOCL_LAPACK_FUNC( dtgexc , DTGEXC )
#define F77_dtgsen AOCL_LAPACK_FUNC( dtgsen , DTGSEN )
#define F77_dtgsja AOCL_LAPACK_FUNC( dtgsja , DTGSJA )
#define F77_dtgsna AOCL_LAPACK_FUNC( dtgsna , DTGSNA )
#define F77_dtgsy2 AOCL_LAPACK_FUNC( dtgsy2 , DTGSY2 )
#define F77_dtgsyl AOCL_LAPACK_FUNC( dtgsyl , DTGSYL )
#define F77_dtpcon AOCL_LAPACK_FUNC( dtpcon , DTPCON )
#define F77_dtplqt2 AOCL_LAPACK_FUNC( dtplqt2 , DTPLQT2 )
#define F77_dtplqt AOCL_LAPACK_FUNC( dtplqt , DTPLQT )
#define F77_dtpmlqt AOCL_LAPACK_FUNC( dtpmlqt , DTPMLQT )
#define F77_dtpmqrt AOCL_LAPACK_FUNC( dtpmqrt , DTPMQRT )
#define F77_dtpqrt2 AOCL_LAPACK_FUNC( dtpqrt2 , DTPQRT2 )
#define F77_dtpqrt AOCL_LAPACK_FUNC( dtpqrt , DTPQRT )
#define F77_dtprfb AOCL_LAPACK_FUNC( dtprfb , DTPRFB )
#define F77_dtprfs AOCL_LAPACK_FUNC( dtprfs , DTPRFS )
#define F77_dtptri AOCL_LAPACK_FUNC( dtptri , DTPTRI )
#define F77_dtptrs AOCL_LAPACK_FUNC( dtptrs , DTPTRS )
#define F77_dtpttf AOCL_LAPACK_FUNC( dtpttf , DTPTTF )
#define F77_dtpttr AOCL_LAPACK_FUNC( dtpttr , DTPTTR )
#define F77_dtrcon AOCL_LAPACK_FUNC( dtrcon , DTRCON )
#define F77_dtrevc3 AOCL_LAPACK_FUNC( dtrevc3 , DTREVC3 )
#define F77_dtrevc AOCL_LAPACK_FUNC( dtrevc , DTREVC )
#define F77_dtrexc AOCL_LAPACK_FUNC( dtrexc , DTREXC )
#define F77_dtrrfs AOCL_LAPACK_FUNC( dtrrfs , DTRRFS )
#define F77_dtrsen AOCL_LAPACK_FUNC( dtrsen , DTRSEN )
#define F77_dtrsna AOCL_LAPACK_FUNC( dtrsna , DTRSNA )
#define F77_dtrsyl AOCL_LAPACK_FUNC( dtrsyl , DTRSYL )
#define F77_dtrti2 AOCL_LAPACK_FUNC( dtrti2 , DTRTI2 )
#define F77_dtrtri AOCL_LAPACK_FUNC( dtrtri , DTRTRI )
#define F77_dtrtrs AOCL_LAPACK_FUNC( dtrtrs , DTRTRS )
#define F77_dtrttf AOCL_LAPACK_FUNC( dtrttf , DTRTTF )
#define F77_dtrttp AOCL_LAPACK_FUNC( dtrttp , DTRTTP )
#define F77_dtzrqf AOCL_LAPACK_FUNC( dtzrqf , DTZRQF )
#define F77_dtzrzf AOCL_LAPACK_FUNC( dtzrzf , DTZRZF )
#define F77_dzsum1 AOCL_LAPACK_FUNC( dzsum1 , DZSUM1 )
#define F77_icmax1 AOCL_LAPACK_FUNC( icmax1 , ICMAX1 )
#define F77_ilaclc AOCL_LAPACK_FUNC( ilaclc , ILACLC )
#define F77_ilaclr AOCL_LAPACK_FUNC( ilaclr , ILACLR )
#define F77_iladiag F77_FUNC( iladiag , ILADIAG )
#define F77_iladlc AOCL_LAPACK_FUNC( iladlc , ILADLC )
#define F77_iladlr AOCL_LAPACK_FUNC( iladlr , ILADLR )
#define F77_ilaenv2stage AOCL_LAPACK_FUNC( ilaenv2stage , ILAENV2STAGE )
#define F77_ilaprec F77_FUNC( ilaprec , ILAPREC )
#define F77_ilaslc AOCL_LAPACK_FUNC( ilaslc , ILASLC )
#define F77_ilaslr AOCL_LAPACK_FUNC( ilaslr , ILASLR )
#define F77_ilatrans F77_FUNC( ilatrans , ILATRANS )
#define F77_ilauplo F77_FUNC( ilauplo , ILAUPLO )
#define F77_ilaver F77_FUNC( ilaver , ILAVER )
#define F77_ilazlc AOCL_LAPACK_FUNC( ilazlc , ILAZLC )
#define F77_ilazlr AOCL_LAPACK_FUNC( ilazlr , ILAZLR )
#define F77_iparam2stage F77_FUNC( iparam2stage , IPARAM2STAGE )
#define F77_izmax1 AOCL_LAPACK_FUNC( izmax1 , IZMAX1 )
#define F77_sbbcsd AOCL_LAPACK_FUNC( sbbcsd , SBBCSD )
#define F77_sbdsdc AOCL_LAPACK_FUNC( sbdsdc , SBDSDC )
#define F77_sbdsqr AOCL_LAPACK_FUNC( sbdsqr , SBDSQR )
#define F77_sbdsvdx AOCL_LAPACK_FUNC( sbdsvdx , SBDSVDX )
#define F77_scombssq F77_FUNC( scombssq , SCOMBSSQ )
#define F77_scsum1 AOCL_LAPACK_FUNC( scsum1 , SCSUM1 )
#define F77_sdisna AOCL_LAPACK_FUNC( sdisna , SDISNA )
#define F77_sgbbrd AOCL_LAPACK_FUNC( sgbbrd , SGBBRD )
#define F77_sgbcon AOCL_LAPACK_FUNC( sgbcon , SGBCON )
#define F77_sgbequb AOCL_LAPACK_FUNC( sgbequb , SGBEQUB )
#define F77_sgbequ AOCL_LAPACK_FUNC( sgbequ , SGBEQU )
#define F77_sgbrfs AOCL_LAPACK_FUNC( sgbrfs , SGBRFS )
#define F77_sgbrfsx F77_FUNC( sgbrfsx , SGBRFSX )
#define F77_sgbsv AOCL_LAPACK_FUNC( sgbsv , SGBSV )
#define F77_sgbsvx AOCL_LAPACK_FUNC( sgbsvx , SGBSVX )
#define F77_sgbsvxx F77_FUNC( sgbsvxx , SGBSVXX )
#define F77_sgbtf2 AOCL_LAPACK_FUNC( sgbtf2 , SGBTF2 )
#define F77_sgbtrf AOCL_LAPACK_FUNC( sgbtrf , SGBTRF )
#define F77_sgbtrs AOCL_LAPACK_FUNC( sgbtrs , SGBTRS )
#define F77_sgebak AOCL_LAPACK_FUNC( sgebak , SGEBAK )
#define F77_sgebal AOCL_LAPACK_FUNC( sgebal , SGEBAL )
#define F77_sgebd2 AOCL_LAPACK_FUNC( sgebd2 , SGEBD2 )
#define F77_sgebrd AOCL_LAPACK_FUNC( sgebrd , SGEBRD )
#define F77_sgecon AOCL_LAPACK_FUNC( sgecon , SGECON )
#define F77_sgeequb AOCL_LAPACK_FUNC( sgeequb , SGEEQUB )
#define F77_sgeequ AOCL_LAPACK_FUNC( sgeequ , SGEEQU )
#define F77_sgees AOCL_LAPACK_FUNC( sgees , SGEES )
#define F77_sgeesx AOCL_LAPACK_FUNC( sgeesx , SGEESX )
#define F77_sgeev AOCL_LAPACK_FUNC( sgeev , SGEEV )
#define F77_sgeevx AOCL_LAPACK_FUNC( sgeevx , SGEEVX )
#define F77_sgegs AOCL_LAPACK_FUNC( sgegs , SGEGS )
#define F77_sgegv AOCL_LAPACK_FUNC( sgegv , SGEGV )
#define F77_sgehd2 AOCL_LAPACK_FUNC( sgehd2 , SGEHD2 )
#define F77_sgehrd AOCL_LAPACK_FUNC( sgehrd , SGEHRD )
#define F77_sgejsv AOCL_LAPACK_FUNC( sgejsv , SGEJSV )
#define F77_sgelq2 AOCL_LAPACK_FUNC( sgelq2 , SGELQ2 )
#define F77_sgelq AOCL_LAPACK_FUNC( sgelq , SGELQ )
#define F77_sgelqf AOCL_LAPACK_FUNC( sgelqf , SGELQF )
#define F77_sgelqt3 AOCL_LAPACK_FUNC( sgelqt3 , SGELQT3 )
#define F77_sgelqt AOCL_LAPACK_FUNC( sgelqt , SGELQT )
#define F77_sgelsd AOCL_LAPACK_FUNC( sgelsd , SGELSD )
#define F77_sgels AOCL_LAPACK_FUNC( sgels , SGELS )
#define F77_sgelss AOCL_LAPACK_FUNC( sgelss , SGELSS )
#define F77_sgelsx AOCL_LAPACK_FUNC( sgelsx , SGELSX )
#define F77_sgelsy AOCL_LAPACK_FUNC( sgelsy , SGELSY )
#define F77_sgemlq AOCL_LAPACK_FUNC( sgemlq , SGEMLQ )
#define F77_sgemlqt AOCL_LAPACK_FUNC( sgemlqt , SGEMLQT )
#define F77_sgemqr AOCL_LAPACK_FUNC( sgemqr , SGEMQR )
#define F77_sgemqrt AOCL_LAPACK_FUNC( sgemqrt , SGEMQRT )
#define F77_sgeql2 AOCL_LAPACK_FUNC( sgeql2 , SGEQL2 )
#define F77_sgeqlf AOCL_LAPACK_FUNC( sgeqlf , SGEQLF )
#define F77_sgeqp3 AOCL_LAPACK_FUNC( sgeqp3 , SGEQP3 )
#define F77_sgeqpf AOCL_LAPACK_FUNC( sgeqpf , SGEQPF )
#define F77_sgeqr2 AOCL_LAPACK_FUNC( sgeqr2 , SGEQR2 )
#define F77_sgeqr2p AOCL_LAPACK_FUNC( sgeqr2p , SGEQR2P )
#define F77_sgeqr AOCL_LAPACK_FUNC( sgeqr , SGEQR )
#define F77_sgeqrf AOCL_LAPACK_FUNC( sgeqrf , SGEQRF )
#define F77_sgeqrfp AOCL_LAPACK_FUNC( sgeqrfp , SGEQRFP )
#define F77_sgeqrt2 AOCL_LAPACK_FUNC( sgeqrt2 , SGEQRT2 )
#define F77_sgeqrt3 AOCL_LAPACK_FUNC( sgeqrt3 , SGEQRT3 )
#define F77_sgeqrt AOCL_LAPACK_FUNC( sgeqrt , SGEQRT )
#define F77_sgerfs AOCL_LAPACK_FUNC( sgerfs , SGERFS )
#define F77_sgerfsx F77_FUNC( sgerfsx , SGERFSX )
#define F77_sgerq2 AOCL_LAPACK_FUNC( sgerq2 , SGERQ2 )
#define F77_sgerqf AOCL_LAPACK_FUNC( sgerqf , SGERQF )
#define F77_sgesc2 AOCL_LAPACK_FUNC( sgesc2 , SGESC2 )
#define F77_sgesdd AOCL_LAPACK_FUNC( sgesdd , SGESDD )
#define F77_sgesvd AOCL_LAPACK_FUNC( sgesvd , SGESVD )
#define F77_sgesvdq AOCL_LAPACK_FUNC( sgesvdq , SGESVDQ )
#define F77_sgesvdx AOCL_LAPACK_FUNC( sgesvdx , SGESVDX )
#define F77_sgesv AOCL_LAPACK_FUNC( sgesv , SGESV )
#define F77_sgesvj AOCL_LAPACK_FUNC( sgesvj , SGESVJ )
#define F77_sgesvx AOCL_LAPACK_FUNC( sgesvx , SGESVX )
#define F77_sgesvxx F77_FUNC( sgesvxx , SGESVXX )
#define F77_sgetc2 AOCL_LAPACK_FUNC( sgetc2 , SGETC2 )
#define F77_sgetf2 AOCL_LAPACK_FUNC( sgetf2 , SGETF2 )
#define F77_sgetrf2 AOCL_LAPACK_FUNC( sgetrf2 , SGETRF2 )
#define F77_sgetrf AOCL_LAPACK_FUNC( sgetrf , SGETRF )
#define F77_sgetri AOCL_LAPACK_FUNC( sgetri , SGETRI )
#define F77_sgetrs AOCL_LAPACK_FUNC( sgetrs , SGETRS )
#define F77_sgetsls AOCL_LAPACK_FUNC( sgetsls , SGETSLS )
#define F77_sggbak AOCL_LAPACK_FUNC( sggbak , SGGBAK )
#define F77_sggbal AOCL_LAPACK_FUNC( sggbal , SGGBAL )
#define F77_sgges3 AOCL_LAPACK_FUNC( sgges3 , SGGES3 )
#define F77_sgges AOCL_LAPACK_FUNC( sgges , SGGES )
#define F77_sggesx AOCL_LAPACK_FUNC( sggesx , SGGESX )
#define F77_sggev3 AOCL_LAPACK_FUNC( sggev3 , SGGEV3 )
#define F77_sggev AOCL_LAPACK_FUNC( sggev , SGGEV )
#define F77_sggevx AOCL_LAPACK_FUNC( sggevx , SGGEVX )
#define F77_sggglm AOCL_LAPACK_FUNC( sggglm , SGGGLM )
#define F77_sgghd3 AOCL_LAPACK_FUNC( sgghd3 , SGGHD3 )
#define F77_sgghrd AOCL_LAPACK_FUNC( sgghrd , SGGHRD )
#define F77_sgglse AOCL_LAPACK_FUNC( sgglse , SGGLSE )
#define F77_sggqrf AOCL_LAPACK_FUNC( sggqrf , SGGQRF )
#define F77_sggrqf AOCL_LAPACK_FUNC( sggrqf , SGGRQF )
#define F77_sggsvd3 AOCL_LAPACK_FUNC( sggsvd3 , SGGSVD3 )
#define F77_sggsvd AOCL_LAPACK_FUNC( sggsvd , SGGSVD )
#define F77_sggsvp3 AOCL_LAPACK_FUNC( sggsvp3 , SGGSVP3 )
#define F77_sggsvp AOCL_LAPACK_FUNC( sggsvp , SGGSVP )
#define F77_sgsvj0 AOCL_LAPACK_FUNC( sgsvj0 , SGSVJ0 )
#define F77_sgsvj1 AOCL_LAPACK_FUNC( sgsvj1 , SGSVJ1 )
#define F77_sgtcon AOCL_LAPACK_FUNC( sgtcon , SGTCON )
#define F77_sgtrfs AOCL_LAPACK_FUNC( sgtrfs , SGTRFS )
#define F77_sgtsv AOCL_LAPACK_FUNC( sgtsv , SGTSV )
#define F77_sgtsvx AOCL_LAPACK_FUNC( sgtsvx , SGTSVX )
#define F77_sgttrf AOCL_LAPACK_FUNC( sgttrf , SGTTRF )
#define F77_sgttrs AOCL_LAPACK_FUNC( sgttrs , SGTTRS )
#define F77_sgtts2 AOCL_LAPACK_FUNC( sgtts2 , SGTTS2 )
#define F77_shgeqz AOCL_LAPACK_FUNC( shgeqz , SHGEQZ )
#define F77_shsein AOCL_LAPACK_FUNC( shsein , SHSEIN )
#define F77_shseqr AOCL_LAPACK_FUNC( shseqr , SHSEQR )
#define F77_sisnan F77_FUNC( sisnan , SISNAN )
#define F77_slabad F77_FUNC( slabad , SLABAD )
#define F77_slabrd AOCL_LAPACK_FUNC( slabrd , SLABRD )
#define F77_slacn2 AOCL_LAPACK_FUNC( slacn2 , SLACN2 )
#define F77_slacon AOCL_LAPACK_FUNC( slacon , SLACON )
#define F77_slacpy AOCL_LAPACK_FUNC( slacpy , SLACPY )
#define F77_sladiv F77_FUNC( sladiv , SLADIV )
#define F77_slae2 F77_FUNC( slae2 , SLAE2 )
#define F77_slaebz AOCL_LAPACK_FUNC( slaebz , SLAEBZ )
#define F77_slaed0 AOCL_LAPACK_FUNC( slaed0 , SLAED0 )
#define F77_slaed1 AOCL_LAPACK_FUNC( slaed1 , SLAED1 )
#define F77_slaed2 AOCL_LAPACK_FUNC( slaed2 , SLAED2 )
#define F77_slaed3 AOCL_LAPACK_FUNC( slaed3 , SLAED3 )
#define F77_slaed4 AOCL_LAPACK_FUNC( slaed4 , SLAED4 )
#define F77_slaed5 AOCL_LAPACK_FUNC( slaed5 , SLAED5 )
#define F77_slaed6 AOCL_LAPACK_FUNC( slaed6 , SLAED6 )
#define F77_slaed7 AOCL_LAPACK_FUNC( slaed7 , SLAED7 )
#define F77_slaed8 AOCL_LAPACK_FUNC( slaed8 , SLAED8 )
#define F77_slaed9 AOCL_LAPACK_FUNC( slaed9 , SLAED9 )
#define F77_slaeda AOCL_LAPACK_FUNC( slaeda , SLAEDA )
#define F77_slaein AOCL_LAPACK_FUNC( slaein , SLAEIN )
#define F77_slaev2 F77_FUNC( slaev2 , SLAEV2 )
#define F77_slaexc AOCL_LAPACK_FUNC( slaexc , SLAEXC )
#define F77_slag2d AOCL_LAPACK_FUNC( slag2d , SLAG2D )
#define F77_slag2 AOCL_LAPACK_FUNC( slag2 , SLAG2 )
#define F77_sla_gbamv F77_FUNC( sla_gbamv , SLA_GBAMV )
#define F77_sla_gbrcond F77_FUNC( sla_gbrcond , SLA_GBRCOND )
#define F77_sla_gbrfsx_extended F77_FUNC( sla_gbrfsx_extended , SLA_GBRFSX_EXTENDED )
#define F77_sla_gbrpvgrw F77_FUNC( sla_gbrpvgrw , SLA_GBRPVGRW )
#define F77_sla_geamv F77_FUNC( sla_geamv , SLA_GEAMV )
#define F77_sla_gercond F77_FUNC( sla_gercond , SLA_GERCOND )
#define F77_sla_gerfsx_extended F77_FUNC( sla_gerfsx_extended , SLA_GERFSX_EXTENDED )
#define F77_sla_gerpvgrw F77_FUNC( sla_gerpvgrw , SLA_GERPVGRW )
#define F77_slags2 F77_FUNC( slags2 , SLAGS2 )
#define F77_slagtf AOCL_LAPACK_FUNC( slagtf , SLAGTF )
#define F77_slagtm AOCL_LAPACK_FUNC( slagtm , SLAGTM )
#define F77_slagts AOCL_LAPACK_FUNC( slagts , SLAGTS )
#define F77_slagv2 AOCL_LAPACK_FUNC( slagv2 , SLAGV2 )
#define F77_slahqr AOCL_LAPACK_FUNC( slahqr , SLAHQR )
#define F77_slahr2 AOCL_LAPACK_FUNC( slahr2 , SLAHR2 )
#define F77_slahrd AOCL_LAPACK_FUNC( slahrd , SLAHRD )
#define F77_slaic1 AOCL_LAPACK_FUNC( slaic1 , SLAIC1 )
#define F77_slaisnan F77_FUNC( slaisnan , SLAISNAN )
#define F77_sla_lin_berr F77_FUNC( sla_lin_berr , SLA_LIN_BERR )
#define F77_slaln2 AOCL_LAPACK_FUNC( slaln2 , SLALN2 )
#define F77_slals0 AOCL_LAPACK_FUNC( slals0 , SLALS0 )
#define F77_slalsa AOCL_LAPACK_FUNC( slalsa , SLALSA )
#define F77_slalsd AOCL_LAPACK_FUNC( slalsd , SLALSD )
#define F77_slamrg AOCL_LAPACK_FUNC( slamrg , SLAMRG )
#define F77_slamswlq AOCL_LAPACK_FUNC( slamswlq , SLAMSWLQ )
#define F77_slamtsqr AOCL_LAPACK_FUNC( slamtsqr , SLAMTSQR )
#define F77_slaneg AOCL_LAPACK_FUNC( slaneg , SLANEG )
#define F77_slangb AOCL_LAPACK_FUNC( slangb , SLANGB )
#define F77_slange AOCL_LAPACK_FUNC( slange , SLANGE )
#define F77_slangt AOCL_LAPACK_FUNC( slangt , SLANGT )
#define F77_slanhs AOCL_LAPACK_FUNC( slanhs , SLANHS )
#define F77_slansb AOCL_LAPACK_FUNC( slansb , SLANSB )
#define F77_slansf AOCL_LAPACK_FUNC( slansf , SLANSF )
#define F77_slansp AOCL_LAPACK_FUNC( slansp , SLANSP )
#define F77_slanst AOCL_LAPACK_FUNC( slanst , SLANST )
#define F77_slansy AOCL_LAPACK_FUNC( slansy , SLANSY )
#define F77_slantb AOCL_LAPACK_FUNC( slantb , SLANTB )
#define F77_slantp AOCL_LAPACK_FUNC( slantp , SLANTP )
#define F77_slantr AOCL_LAPACK_FUNC( slantr , SLANTR )
#define F77_slanv2 F77_FUNC( slanv2 , SLANV2 )
#define F77_slaorhr_col_getrfnp2 F77_FUNC( slaorhr_col_getrfnp2 , SLAORHR_COL_GETRFNP2 )
#define F77_slaorhr_col_getrfnp F77_FUNC( slaorhr_col_getrfnp , SLAORHR_COL_GETRFNP )
#define F77_slapll AOCL_LAPACK_FUNC( slapll , SLAPLL )
#define F77_slapmr AOCL_LAPACK_FUNC( slapmr , SLAPMR )
#define F77_slapmt AOCL_LAPACK_FUNC( slapmt , SLAPMT )
#define F77_sla_porcond F77_FUNC( sla_porcond , SLA_PORCOND )
#define F77_sla_porfsx_extended F77_FUNC( sla_porfsx_extended , SLA_PORFSX_EXTENDED )
#define F77_sla_porpvgrw F77_FUNC( sla_porpvgrw , SLA_PORPVGRW )
#define F77_slapy2 F77_FUNC( slapy2 , SLAPY2 )
#define F77_slapy3 F77_FUNC( slapy3 , SLAPY3 )
#define F77_slaqgb AOCL_LAPACK_FUNC( slaqgb , SLAQGB )
#define F77_slaqge AOCL_LAPACK_FUNC( slaqge , SLAQGE )
#define F77_slaqp2 AOCL_LAPACK_FUNC( slaqp2 , SLAQP2 )
#define F77_slaqps AOCL_LAPACK_FUNC( slaqps , SLAQPS )
#define F77_slaqr0 AOCL_LAPACK_FUNC( slaqr0 , SLAQR0 )
#define F77_slaqr1 AOCL_LAPACK_FUNC( slaqr1 , SLAQR1 )
#define F77_slaqr2 AOCL_LAPACK_FUNC( slaqr2 , SLAQR2 )
#define F77_slaqr3 AOCL_LAPACK_FUNC( slaqr3 , SLAQR3 )
#define F77_slaqr4 AOCL_LAPACK_FUNC( slaqr4 , SLAQR4 )
#define F77_slaqr5 AOCL_LAPACK_FUNC( slaqr5 , SLAQR5 )
#define F77_slaqsb AOCL_LAPACK_FUNC( slaqsb , SLAQSB )
#define F77_slaqsp AOCL_LAPACK_FUNC( slaqsp , SLAQSP )
#define F77_slaqsy AOCL_LAPACK_FUNC( slaqsy , SLAQSY )
#define F77_slaqtr AOCL_LAPACK_FUNC( slaqtr , SLAQTR )
#define F77_slar1v AOCL_LAPACK_FUNC( slar1v , SLAR1V )
#define F77_slar2v AOCL_LAPACK_FUNC( slar2v , SLAR2V )
#define F77_slarfb AOCL_LAPACK_FUNC( slarfb , SLARFB )
#define F77_slarf AOCL_LAPACK_FUNC( slarf , SLARF )
#define F77_slarfg AOCL_LAPACK_FUNC( slarfg , SLARFG )
#define F77_slarfgp AOCL_LAPACK_FUNC( slarfgp , SLARFGP )
#define F77_slarft AOCL_LAPACK_FUNC( slarft , SLARFT )
#define F77_slarfx AOCL_LAPACK_FUNC( slarfx , SLARFX )
#define F77_slarfy AOCL_LAPACK_FUNC( slarfy , SLARFY )
#define F77_slargv AOCL_LAPACK_FUNC( slargv , SLARGV )
#define F77_slarnv AOCL_LAPACK_FUNC( slarnv , SLARNV )
#define F77_slarra AOCL_LAPACK_FUNC( slarra , SLARRA )
#define F77_slarrb AOCL_LAPACK_FUNC( slarrb , SLARRB )
#define F77_slarrc AOCL_LAPACK_FUNC( slarrc , SLARRC )
#define F77_slarrd AOCL_LAPACK_FUNC( slarrd , SLARRD )
#define F77_slarre AOCL_LAPACK_FUNC( slarre , SLARRE )
#define F77_slarrf AOCL_LAPACK_FUNC( slarrf , SLARRF )
#define F77_slarrj AOCL_LAPACK_FUNC( slarrj , SLARRJ )
#define F77_slarrk AOCL_LAPACK_FUNC( slarrk , SLARRK )
#define F77_slarrr AOCL_LAPACK_FUNC( slarrr , SLARRR )
#define F77_slarrv AOCL_LAPACK_FUNC( slarrv , SLARRV )
#define F77_slarscl2 AOCL_LAPACK_FUNC( slarscl2 , SLARSCL2 )
#define F77_slartg F77_FUNC( slartg , SLARTG )
#define F77_slartgp F77_FUNC( slartgp , SLARTGP )
#define F77_slartgs F77_FUNC( slartgs , SLARTGS )
#define F77_slartv AOCL_LAPACK_FUNC( slartv , SLARTV )
#define F77_slaruv AOCL_LAPACK_FUNC( slaruv , SLARUV )
#define F77_slarzb AOCL_LAPACK_FUNC( slarzb , SLARZB )
#define F77_slarz AOCL_LAPACK_FUNC( slarz , SLARZ )
#define F77_slarzt AOCL_LAPACK_FUNC( slarzt , SLARZT )
#define F77_slas2 F77_FUNC( slas2 , SLAS2 )
#define F77_slascl2 AOCL_LAPACK_FUNC( slascl2 , SLASCL2 )
#define F77_slascl AOCL_LAPACK_FUNC( slascl , SLASCL )
#define F77_slasd0 AOCL_LAPACK_FUNC( slasd0 , SLASD0 )
#define F77_slasd1 AOCL_LAPACK_FUNC( slasd1 , SLASD1 )
#define F77_slasd2 AOCL_LAPACK_FUNC( slasd2 , SLASD2 )
#define F77_slasd3 AOCL_LAPACK_FUNC( slasd3 , SLASD3 )
#define F77_slasd4 AOCL_LAPACK_FUNC( slasd4 , SLASD4 )
#define F77_slasd5 AOCL_LAPACK_FUNC( slasd5 , SLASD5 )
#define F77_slasd6 AOCL_LAPACK_FUNC( slasd6 , SLASD6 )
#define F77_slasd7 AOCL_LAPACK_FUNC( slasd7 , SLASD7 )
#define F77_slasd8 AOCL_LAPACK_FUNC( slasd8 , SLASD8 )
#define F77_slasda AOCL_LAPACK_FUNC( slasda , SLASDA )
#define F77_slasdq AOCL_LAPACK_FUNC( slasdq , SLASDQ )
#define F77_slasdt AOCL_LAPACK_FUNC( slasdt , SLASDT )
#define F77_slaset AOCL_LAPACK_FUNC( slaset , SLASET )
#define F77_slasq1 AOCL_LAPACK_FUNC( slasq1 , SLASQ1 )
#define F77_slasq2 AOCL_LAPACK_FUNC( slasq2 , SLASQ2 )
#define F77_slasq3 AOCL_LAPACK_FUNC( slasq3 , SLASQ3 )
#define F77_slasq4 AOCL_LAPACK_FUNC( slasq4 , SLASQ4 )
#define F77_slasq5 AOCL_LAPACK_FUNC( slasq5 , SLASQ5 )
#define F77_slasq6 AOCL_LAPACK_FUNC( slasq6 , SLASQ6 )
#define F77_slasr AOCL_LAPACK_FUNC( slasr , SLASR )
#define F77_slasrt AOCL_LAPACK_FUNC( slasrt , SLASRT )
#define F77_slassq AOCL_LAPACK_FUNC( slassq , SLASSQ )
#define F77_slasv2 F77_FUNC( slasv2 , SLASV2 )
#define F77_slaswlq AOCL_LAPACK_FUNC( slaswlq , SLASWLQ )
#define F77_slaswp AOCL_LAPACK_FUNC( slaswp , SLASWP )
#define F77_slasy2 AOCL_LAPACK_FUNC( slasy2 , SLASY2 )
#define F77_sla_syamv F77_FUNC( sla_syamv , SLA_SYAMV )
#define F77_slasyf_aa AOCL_LAPACK_FUNC( slasyf_aa , SLASYF_AA )
#define F77_slasyf AOCL_LAPACK_FUNC( slasyf , SLASYF )
#define F77_slasyf_rk AOCL_LAPACK_FUNC( slasyf_rk , SLASYF_RK )
#define F77_slasyf_rook AOCL_LAPACK_FUNC( slasyf_rook , SLASYF_ROOK )
#define F77_sla_syrcond F77_FUNC( sla_syrcond , SLA_SYRCOND )
#define F77_sla_syrfsx_extended F77_FUNC( sla_syrfsx_extended , SLA_SYRFSX_EXTENDED )
#define F77_sla_syrpvgrw F77_FUNC( sla_syrpvgrw , SLA_SYRPVGRW )
#define F77_slatbs AOCL_LAPACK_FUNC( slatbs , SLATBS )
#define F77_slatdf AOCL_LAPACK_FUNC( slatdf , SLATDF )
#define F77_slatps AOCL_LAPACK_FUNC( slatps , SLATPS )
#define F77_slatrd AOCL_LAPACK_FUNC( slatrd , SLATRD )
#define F77_slatrs AOCL_LAPACK_FUNC( slatrs , SLATRS )
#define F77_slatrz AOCL_LAPACK_FUNC( slatrz , SLATRZ )
#define F77_slatsqr AOCL_LAPACK_FUNC( slatsqr , SLATSQR )
#define F77_slatzm AOCL_LAPACK_FUNC( slatzm , SLATZM )
#define F77_slauu2 AOCL_LAPACK_FUNC( slauu2 , SLAUU2 )
#define F77_slauum AOCL_LAPACK_FUNC( slauum , SLAUUM )
#define F77_sla_wwaddw F77_FUNC( sla_wwaddw , SLA_WWADDW )
#define F77_sopmtr AOCL_LAPACK_FUNC( sopmtr , SOPMTR )
#define F77_sorbdb1 AOCL_LAPACK_FUNC( sorbdb1 , SORBDB1 )
#define F77_sorbdb2 AOCL_LAPACK_FUNC( sorbdb2 , SORBDB2 )
#define F77_sorbdb3 AOCL_LAPACK_FUNC( sorbdb3 , SORBDB3 )
#define F77_sorbdb4 AOCL_LAPACK_FUNC( sorbdb4 , SORBDB4 )
#define F77_sorbdb5 AOCL_LAPACK_FUNC( sorbdb5 , SORBDB5 )
#define F77_sorbdb6 AOCL_LAPACK_FUNC( sorbdb6 , SORBDB6 )
#define F77_sorbdb AOCL_LAPACK_FUNC( sorbdb , SORBDB )
#define F77_sorg2l AOCL_LAPACK_FUNC( sorg2l , SORG2L )
#define F77_sorg2r AOCL_LAPACK_FUNC( sorg2r , SORG2R )
#define F77_sorgbr AOCL_LAPACK_FUNC( sorgbr , SORGBR )
#define F77_sorgl2 AOCL_LAPACK_FUNC( sorgl2 , SORGL2 )
#define F77_sorglq AOCL_LAPACK_FUNC( sorglq , SORGLQ )
#define F77_sorgql AOCL_LAPACK_FUNC( sorgql , SORGQL )
#define F77_sorgqr AOCL_LAPACK_FUNC( sorgqr , SORGQR )
#define F77_sorgr2 AOCL_LAPACK_FUNC( sorgr2 , SORGR2 )
#define F77_sorgrq AOCL_LAPACK_FUNC( sorgrq , SORGRQ )
#define F77_sorgtr AOCL_LAPACK_FUNC( sorgtr , SORGTR )
#define F77_sorgtsqr AOCL_LAPACK_FUNC( sorgtsqr , SORGTSQR )
#define F77_sorhr_col F77_FUNC( sorhr_col , SORHR_COL )
#define F77_sorm22 AOCL_LAPACK_FUNC( sorm22 , SORM22 )
#define F77_sorm2l AOCL_LAPACK_FUNC( sorm2l , SORM2L )
#define F77_sorm2r AOCL_LAPACK_FUNC( sorm2r , SORM2R )
#define F77_sormbr AOCL_LAPACK_FUNC( sormbr , SORMBR )
#define F77_sorml2 AOCL_LAPACK_FUNC( sorml2 , SORML2 )
#define F77_sormlq AOCL_LAPACK_FUNC( sormlq , SORMLQ )
#define F77_sormql AOCL_LAPACK_FUNC( sormql , SORMQL )
#define F77_sormqr AOCL_LAPACK_FUNC( sormqr , SORMQR )
#define F77_sormr2 AOCL_LAPACK_FUNC( sormr2 , SORMR2 )
#define F77_sormr3 AOCL_LAPACK_FUNC( sormr3 , SORMR3 )
#define F77_sormrq AOCL_LAPACK_FUNC( sormrq , SORMRQ )
#define F77_sormrz AOCL_LAPACK_FUNC( sormrz , SORMRZ )
#define F77_sormtr AOCL_LAPACK_FUNC( sormtr , SORMTR )
#define F77_spbcon AOCL_LAPACK_FUNC( spbcon , SPBCON )
#define F77_spbequ AOCL_LAPACK_FUNC( spbequ , SPBEQU )
#define F77_spbrfs AOCL_LAPACK_FUNC( spbrfs , SPBRFS )
#define F77_spbstf AOCL_LAPACK_FUNC( spbstf , SPBSTF )
#define F77_spbsv AOCL_LAPACK_FUNC( spbsv , SPBSV )
#define F77_spbsvx AOCL_LAPACK_FUNC( spbsvx , SPBSVX )
#define F77_spbtf2 AOCL_LAPACK_FUNC( spbtf2 , SPBTF2 )
#define F77_spbtrf AOCL_LAPACK_FUNC( spbtrf , SPBTRF )
#define F77_spbtrs AOCL_LAPACK_FUNC( spbtrs , SPBTRS )
#define F77_spftrf AOCL_LAPACK_FUNC( spftrf , SPFTRF )
#define F77_spftri AOCL_LAPACK_FUNC( spftri , SPFTRI )
#define F77_spftrs AOCL_LAPACK_FUNC( spftrs , SPFTRS )
#define F77_spocon AOCL_LAPACK_FUNC( spocon , SPOCON )
#define F77_spoequb AOCL_LAPACK_FUNC( spoequb , SPOEQUB )
#define F77_spoequ AOCL_LAPACK_FUNC( spoequ , SPOEQU )
#define F77_sporfs AOCL_LAPACK_FUNC( sporfs , SPORFS )
#define F77_sporfsx F77_FUNC( sporfsx , SPORFSX )
#define F77_sposv AOCL_LAPACK_FUNC( sposv , SPOSV )
#define F77_sposvx AOCL_LAPACK_FUNC( sposvx , SPOSVX )
#define F77_sposvxx F77_FUNC( sposvxx , SPOSVXX )
#define F77_spotf2 AOCL_LAPACK_FUNC( spotf2 , SPOTF2 )
#define F77_spotrf2 AOCL_LAPACK_FUNC( spotrf2 , SPOTRF2 )
#define F77_spotrf AOCL_LAPACK_FUNC( spotrf , SPOTRF )
#define F77_spotri AOCL_LAPACK_FUNC( spotri , SPOTRI )
#define F77_spotrs AOCL_LAPACK_FUNC( spotrs , SPOTRS )
#define F77_sppcon AOCL_LAPACK_FUNC( sppcon , SPPCON )
#define F77_sppequ AOCL_LAPACK_FUNC( sppequ , SPPEQU )
#define F77_spprfs AOCL_LAPACK_FUNC( spprfs , SPPRFS )
#define F77_sppsv AOCL_LAPACK_FUNC( sppsv , SPPSV )
#define F77_sppsvx AOCL_LAPACK_FUNC( sppsvx , SPPSVX )
#define F77_spptrf AOCL_LAPACK_FUNC( spptrf , SPPTRF )
#define F77_spptri AOCL_LAPACK_FUNC( spptri , SPPTRI )
#define F77_spptrs AOCL_LAPACK_FUNC( spptrs , SPPTRS )
#define F77_spstf2 AOCL_LAPACK_FUNC( spstf2 , SPSTF2 )
#define F77_spstrf AOCL_LAPACK_FUNC( spstrf , SPSTRF )
#define F77_sptcon AOCL_LAPACK_FUNC( sptcon , SPTCON )
#define F77_spteqr AOCL_LAPACK_FUNC( spteqr , SPTEQR )
#define F77_sptrfs AOCL_LAPACK_FUNC( sptrfs , SPTRFS )
#define F77_sptsv AOCL_LAPACK_FUNC( sptsv , SPTSV )
#define F77_sptsvx AOCL_LAPACK_FUNC( sptsvx , SPTSVX )
#define F77_spttrf AOCL_LAPACK_FUNC( spttrf , SPTTRF )
#define F77_spttrs AOCL_LAPACK_FUNC( spttrs , SPTTRS )
#define F77_sptts2 AOCL_LAPACK_FUNC( sptts2 , SPTTS2 )
#define F77_srscl AOCL_LAPACK_FUNC( srscl , SRSCL )
#define F77_ssb2st_kernels F77_FUNC( ssb2st_kernels , SSB2ST_KERNELS )
#define F77_ssbev_2stage AOCL_LAPACK_FUNC( ssbev_2stage , SSBEV_2STAGE )
#define F77_ssbevd_2stage AOCL_LAPACK_FUNC( ssbevd_2stage , SSBEVD_2STAGE )
#define F77_ssbevd AOCL_LAPACK_FUNC( ssbevd , SSBEVD )
#define F77_ssbev AOCL_LAPACK_FUNC( ssbev , SSBEV )
#define F77_ssbevx_2stage AOCL_LAPACK_FUNC( ssbevx_2stage , SSBEVX_2STAGE )
#define F77_ssbevx AOCL_LAPACK_FUNC( ssbevx , SSBEVX )
#define F77_ssbgst AOCL_LAPACK_FUNC( ssbgst , SSBGST )
#define F77_ssbgvd AOCL_LAPACK_FUNC( ssbgvd , SSBGVD )
#define F77_ssbgv AOCL_LAPACK_FUNC( ssbgv , SSBGV )
#define F77_ssbgvx AOCL_LAPACK_FUNC( ssbgvx , SSBGVX )
#define F77_ssbtrd AOCL_LAPACK_FUNC( ssbtrd , SSBTRD )
#define F77_ssfrk AOCL_LAPACK_FUNC( ssfrk , SSFRK )
#define F77_sspcon AOCL_LAPACK_FUNC( sspcon , SSPCON )
#define F77_sspevd AOCL_LAPACK_FUNC( sspevd , SSPEVD )
#define F77_sspev AOCL_LAPACK_FUNC( sspev , SSPEV )
#define F77_sspevx AOCL_LAPACK_FUNC( sspevx , SSPEVX )
#define F77_sspgst AOCL_LAPACK_FUNC( sspgst , SSPGST )
#define F77_sspgvd AOCL_LAPACK_FUNC( sspgvd , SSPGVD )
#define F77_sspgv AOCL_LAPACK_FUNC( sspgv , SSPGV )
#define F77_sspgvx AOCL_LAPACK_FUNC( sspgvx , SSPGVX )
#define F77_ssprfs AOCL_LAPACK_FUNC( ssprfs , SSPRFS )
#define F77_sspsv AOCL_LAPACK_FUNC( sspsv , SSPSV )
#define F77_sspsvx AOCL_LAPACK_FUNC( sspsvx , SSPSVX )
#define F77_ssptrd AOCL_LAPACK_FUNC( ssptrd , SSPTRD )
#define F77_ssptrf AOCL_LAPACK_FUNC( ssptrf , SSPTRF )
#define F77_ssptri AOCL_LAPACK_FUNC( ssptri , SSPTRI )
#define F77_ssptrs AOCL_LAPACK_FUNC( ssptrs , SSPTRS )
#define F77_sstebz AOCL_LAPACK_FUNC( sstebz , SSTEBZ )
#define F77_sstedc AOCL_LAPACK_FUNC( sstedc , SSTEDC )
#define F77_sstegr AOCL_LAPACK_FUNC( sstegr , SSTEGR )
#define F77_sstein AOCL_LAPACK_FUNC( sstein , SSTEIN )
#define F77_sstemr AOCL_LAPACK_FUNC( sstemr , SSTEMR )
#define F77_ssteqr AOCL_LAPACK_FUNC( ssteqr , SSTEQR )
#define F77_ssterf AOCL_LAPACK_FUNC( ssterf , SSTERF )
#define F77_sstevd AOCL_LAPACK_FUNC( sstevd , SSTEVD )
#define F77_sstev AOCL_LAPACK_FUNC( sstev , SSTEV )
#define F77_sstevr AOCL_LAPACK_FUNC( sstevr , SSTEVR )
#define F77_sstevx AOCL_LAPACK_FUNC( sstevx , SSTEVX )
#define F77_ssycon_3 AOCL_LAPACK_FUNC( ssycon_3 , SSYCON_3 )
#define F77_ssycon AOCL_LAPACK_FUNC( ssycon , SSYCON )
#define F77_ssycon_rook AOCL_LAPACK_FUNC( ssycon_rook , SSYCON_ROOK )
#define F77_ssyconv AOCL_LAPACK_FUNC( ssyconv , SSYCONV )
#define F77_ssyconvf AOCL_LAPACK_FUNC( ssyconvf , SSYCONVF )
#define F77_ssyconvf_rook AOCL_LAPACK_FUNC( ssyconvf_rook , SSYCONVF_ROOK )
#define F77_ssyequb AOCL_LAPACK_FUNC( ssyequb , SSYEQUB )
#define F77_ssyev_2stage AOCL_LAPACK_FUNC( ssyev_2stage , SSYEV_2STAGE )
#define F77_ssyevd_2stage AOCL_LAPACK_FUNC( ssyevd_2stage , SSYEVD_2STAGE )
#define F77_ssyevd AOCL_LAPACK_FUNC( ssyevd , SSYEVD )
#define F77_ssyev  AOCL_LAPACK_FUNC( ssyev  , SSYEV  )
#define F77_ssyev AOCL_LAPACK_FUNC( ssyev , SSYEV )
#define F77_ssyevr_2stage AOCL_LAPACK_FUNC( ssyevr_2stage , SSYEVR_2STAGE )
#define F77_ssyevr AOCL_LAPACK_FUNC( ssyevr , SSYEVR )
#define F77_ssyevx_2stage AOCL_LAPACK_FUNC( ssyevx_2stage , SSYEVX_2STAGE )
#define F77_ssyevx AOCL_LAPACK_FUNC( ssyevx , SSYEVX )
#define F77_ssygs2 AOCL_LAPACK_FUNC( ssygs2 , SSYGS2 )
#define F77_ssygst AOCL_LAPACK_FUNC( ssygst , SSYGST )
#define F77_ssygv_2stage AOCL_LAPACK_FUNC( ssygv_2stage , SSYGV_2STAGE )
#define F77_ssygvd AOCL_LAPACK_FUNC( ssygvd , SSYGVD )
#define F77_ssygv AOCL_LAPACK_FUNC( ssygv , SSYGV )
#define F77_ssygvx AOCL_LAPACK_FUNC( ssygvx , SSYGVX )
#define F77_ssyrfs AOCL_LAPACK_FUNC( ssyrfs , SSYRFS )
#define F77_ssyrfsx F77_FUNC( ssyrfsx , SSYRFSX )
#define F77_ssysv_aa_2stage AOCL_LAPACK_FUNC( ssysv_aa_2stage , SSYSV_AA_2STAGE )
#define F77_ssysv_aa AOCL_LAPACK_FUNC( ssysv_aa , SSYSV_AA )
#define F77_ssysv AOCL_LAPACK_FUNC( ssysv , SSYSV )
#define F77_ssysv_rk AOCL_LAPACK_FUNC( ssysv_rk , SSYSV_RK )
#define F77_ssysv_rook AOCL_LAPACK_FUNC( ssysv_rook , SSYSV_ROOK )
#define F77_ssysvx AOCL_LAPACK_FUNC( ssysvx , SSYSVX )
#define F77_ssysvxx F77_FUNC( ssysvxx , SSYSVXX )
#define F77_ssyswapr AOCL_LAPACK_FUNC( ssyswapr , SSYSWAPR )
#define F77_ssytd2 AOCL_LAPACK_FUNC( ssytd2 , SSYTD2 )
#define F77_ssytf2 AOCL_LAPACK_FUNC( ssytf2 , SSYTF2 )
#define F77_ssytf2_rk AOCL_LAPACK_FUNC( ssytf2_rk , SSYTF2_RK )
#define F77_ssytf2_rook AOCL_LAPACK_FUNC( ssytf2_rook , SSYTF2_ROOK )
#define F77_ssytrd_2stage AOCL_LAPACK_FUNC( ssytrd_2stage , SSYTRD_2STAGE )
#define F77_ssytrd AOCL_LAPACK_FUNC( ssytrd , SSYTRD )
#define F77_ssytrd_sb2st AOCL_LAPACK_FUNC( ssytrd_sb2st , SSYTRD_SB2ST )
#define F77_ssytrd_sy2sb AOCL_LAPACK_FUNC( ssytrd_sy2sb , SSYTRD_SY2SB )
#define F77_ssytrf_aa_2stage AOCL_LAPACK_FUNC( ssytrf_aa_2stage , SSYTRF_AA_2STAGE )
#define F77_ssytrf_aa AOCL_LAPACK_FUNC( ssytrf_aa , SSYTRF_AA )
#define F77_ssytrf AOCL_LAPACK_FUNC( ssytrf , SSYTRF )
#define F77_ssytrf_rk AOCL_LAPACK_FUNC( ssytrf_rk , SSYTRF_RK )
#define F77_ssytrf_rook AOCL_LAPACK_FUNC( ssytrf_rook , SSYTRF_ROOK )
#define F77_ssytri2 AOCL_LAPACK_FUNC( ssytri2 , SSYTRI2 )
#define F77_ssytri2x AOCL_LAPACK_FUNC( ssytri2x , SSYTRI2X )
#define F77_ssytri_3 AOCL_LAPACK_FUNC( ssytri_3 , SSYTRI_3 )
#define F77_ssytri_3x AOCL_LAPACK_FUNC( ssytri_3x , SSYTRI_3X )
#define F77_ssytri AOCL_LAPACK_FUNC( ssytri , SSYTRI )
#define F77_ssytri_rook AOCL_LAPACK_FUNC( ssytri_rook , SSYTRI_ROOK )
#define F77_ssytrs2 AOCL_LAPACK_FUNC( ssytrs2 , SSYTRS2 )
#define F77_ssytrs_3 AOCL_LAPACK_FUNC( ssytrs_3 , SSYTRS_3 )
#define F77_ssytrs_aa_2stage AOCL_LAPACK_FUNC( ssytrs_aa_2stage , SSYTRS_AA_2STAGE )
#define F77_ssytrs_aa AOCL_LAPACK_FUNC( ssytrs_aa , SSYTRS_AA )
#define F77_ssytrs AOCL_LAPACK_FUNC( ssytrs , SSYTRS )
#define F77_ssytrs_rook AOCL_LAPACK_FUNC( ssytrs_rook , SSYTRS_ROOK )
#define F77_stbcon AOCL_LAPACK_FUNC( stbcon , STBCON )
#define F77_stbrfs AOCL_LAPACK_FUNC( stbrfs , STBRFS )
#define F77_stbtrs AOCL_LAPACK_FUNC( stbtrs , STBTRS )
#define F77_stfsm AOCL_LAPACK_FUNC( stfsm , STFSM )
#define F77_stftri AOCL_LAPACK_FUNC( stftri , STFTRI )
#define F77_stfttp AOCL_LAPACK_FUNC( stfttp , STFTTP )
#define F77_stfttr AOCL_LAPACK_FUNC( stfttr , STFTTR )
#define F77_stgevc AOCL_LAPACK_FUNC( stgevc , STGEVC )
#define F77_stgex2 AOCL_LAPACK_FUNC( stgex2 , STGEX2 )
#define F77_stgexc AOCL_LAPACK_FUNC( stgexc , STGEXC )
#define F77_stgsen AOCL_LAPACK_FUNC( stgsen , STGSEN )
#define F77_stgsja AOCL_LAPACK_FUNC( stgsja , STGSJA )
#define F77_stgsna AOCL_LAPACK_FUNC( stgsna , STGSNA )
#define F77_stgsy2 AOCL_LAPACK_FUNC( stgsy2 , STGSY2 )
#define F77_stgsyl AOCL_LAPACK_FUNC( stgsyl , STGSYL )
#define F77_stpcon AOCL_LAPACK_FUNC( stpcon , STPCON )
#define F77_stplqt2 AOCL_LAPACK_FUNC( stplqt2 , STPLQT2 )
#define F77_stplqt AOCL_LAPACK_FUNC( stplqt , STPLQT )
#define F77_stpmlqt AOCL_LAPACK_FUNC( stpmlqt , STPMLQT )
#define F77_stpmqrt AOCL_LAPACK_FUNC( stpmqrt , STPMQRT )
#define F77_stpqrt2 AOCL_LAPACK_FUNC( stpqrt2 , STPQRT2 )
#define F77_stpqrt AOCL_LAPACK_FUNC( stpqrt , STPQRT )
#define F77_stprfb AOCL_LAPACK_FUNC( stprfb , STPRFB )
#define F77_stprfs AOCL_LAPACK_FUNC( stprfs , STPRFS )
#define F77_stptri AOCL_LAPACK_FUNC( stptri , STPTRI )
#define F77_stptrs AOCL_LAPACK_FUNC( stptrs , STPTRS )
#define F77_stpttf AOCL_LAPACK_FUNC( stpttf , STPTTF )
#define F77_stpttr AOCL_LAPACK_FUNC( stpttr , STPTTR )
#define F77_strcon AOCL_LAPACK_FUNC( strcon , STRCON )
#define F77_strevc3 AOCL_LAPACK_FUNC( strevc3 , STREVC3 )
#define F77_strevc AOCL_LAPACK_FUNC( strevc , STREVC )
#define F77_strexc AOCL_LAPACK_FUNC( strexc , STREXC )
#define F77_strrfs AOCL_LAPACK_FUNC( strrfs , STRRFS )
#define F77_strsen AOCL_LAPACK_FUNC( strsen , STRSEN )
#define F77_strsna AOCL_LAPACK_FUNC( strsna , STRSNA )
#define F77_strsyl AOCL_LAPACK_FUNC( strsyl , STRSYL )
#define F77_strti2 AOCL_LAPACK_FUNC( strti2 , STRTI2 )
#define F77_strtri AOCL_LAPACK_FUNC( strtri , STRTRI )
#define F77_strtrs AOCL_LAPACK_FUNC( strtrs , STRTRS )
#define F77_strttf AOCL_LAPACK_FUNC( strttf , STRTTF )
#define F77_strttp AOCL_LAPACK_FUNC( strttp , STRTTP )
#define F77_stzrqf AOCL_LAPACK_FUNC( stzrqf , STZRQF )
#define F77_stzrzf AOCL_LAPACK_FUNC( stzrzf , STZRZF )
#define F77_zbbcsd AOCL_LAPACK_FUNC( zbbcsd , ZBBCSD )
#define F77_zbdsqr AOCL_LAPACK_FUNC( zbdsqr , ZBDSQR )
#define F77_zcgesv AOCL_LAPACK_FUNC( zcgesv , ZCGESV )
#define F77_zcposv AOCL_LAPACK_FUNC( zcposv , ZCPOSV )
#define F77_zdrscl AOCL_LAPACK_FUNC( zdrscl , ZDRSCL )
#define F77_zgbbrd AOCL_LAPACK_FUNC( zgbbrd , ZGBBRD )
#define F77_zgbcon AOCL_LAPACK_FUNC( zgbcon , ZGBCON )
#define F77_zgbequb AOCL_LAPACK_FUNC( zgbequb , ZGBEQUB )
#define F77_zgbequ AOCL_LAPACK_FUNC( zgbequ , ZGBEQU )
#define F77_zgbrfs AOCL_LAPACK_FUNC( zgbrfs , ZGBRFS )
#define F77_zgbrfsx F77_FUNC( zgbrfsx , ZGBRFSX )
#define F77_zgbsv AOCL_LAPACK_FUNC( zgbsv , ZGBSV )
#define F77_zgbsvx AOCL_LAPACK_FUNC( zgbsvx , ZGBSVX )
#define F77_zgbsvxx F77_FUNC( zgbsvxx , ZGBSVXX )
#define F77_zgbtf2 AOCL_LAPACK_FUNC( zgbtf2 , ZGBTF2 )
#define F77_zgbtrf AOCL_LAPACK_FUNC( zgbtrf , ZGBTRF )
#define F77_zgbtrs AOCL_LAPACK_FUNC( zgbtrs , ZGBTRS )
#define F77_zgebak AOCL_LAPACK_FUNC( zgebak , ZGEBAK )
#define F77_zgebal AOCL_LAPACK_FUNC( zgebal , ZGEBAL )
#define F77_zgebd2 AOCL_LAPACK_FUNC( zgebd2 , ZGEBD2 )
#define F77_zgebrd AOCL_LAPACK_FUNC( zgebrd , ZGEBRD )
#define F77_zgecon AOCL_LAPACK_FUNC( zgecon , ZGECON )
#define F77_zgeequb AOCL_LAPACK_FUNC( zgeequb , ZGEEQUB )
#define F77_zgeequ AOCL_LAPACK_FUNC( zgeequ , ZGEEQU )
#define F77_zgees AOCL_LAPACK_FUNC( zgees , ZGEES )
#define F77_zgeesx AOCL_LAPACK_FUNC( zgeesx , ZGEESX )
#define F77_zgeev AOCL_LAPACK_FUNC( zgeev , ZGEEV )
#define F77_zgeevx AOCL_LAPACK_FUNC( zgeevx , ZGEEVX )
#define F77_zgegs AOCL_LAPACK_FUNC( zgegs , ZGEGS )
#define F77_zgegv AOCL_LAPACK_FUNC( zgegv , ZGEGV )
#define F77_zgehd2 AOCL_LAPACK_FUNC( zgehd2 , ZGEHD2 )
#define F77_zgehrd AOCL_LAPACK_FUNC( zgehrd , ZGEHRD )
#define F77_zgejsv AOCL_LAPACK_FUNC( zgejsv , ZGEJSV )
#define F77_zgelq2 AOCL_LAPACK_FUNC( zgelq2 , ZGELQ2 )
#define F77_zgelq AOCL_LAPACK_FUNC( zgelq , ZGELQ )
#define F77_zgelqf AOCL_LAPACK_FUNC( zgelqf , ZGELQF )
#define F77_zgelqt3 AOCL_LAPACK_FUNC( zgelqt3 , ZGELQT3 )
#define F77_zgelqt AOCL_LAPACK_FUNC( zgelqt , ZGELQT )
#define F77_zgelsd AOCL_LAPACK_FUNC( zgelsd , ZGELSD )
#define F77_zgels AOCL_LAPACK_FUNC( zgels , ZGELS )
#define F77_zgelss AOCL_LAPACK_FUNC( zgelss , ZGELSS )
#define F77_zgelsx AOCL_LAPACK_FUNC( zgelsx , ZGELSX )
#define F77_zgelsy AOCL_LAPACK_FUNC( zgelsy , ZGELSY )
#define F77_zgemlq AOCL_LAPACK_FUNC( zgemlq , ZGEMLQ )
#define F77_zgemlqt AOCL_LAPACK_FUNC( zgemlqt , ZGEMLQT )
#define F77_zgemqr AOCL_LAPACK_FUNC( zgemqr , ZGEMQR )
#define F77_zgemqrt AOCL_LAPACK_FUNC( zgemqrt , ZGEMQRT )
#define F77_zgeql2 AOCL_LAPACK_FUNC( zgeql2 , ZGEQL2 )
#define F77_zgeqlf AOCL_LAPACK_FUNC( zgeqlf , ZGEQLF )
#define F77_zgeqp3 AOCL_LAPACK_FUNC( zgeqp3 , ZGEQP3 )
#define F77_zgeqpf AOCL_LAPACK_FUNC( zgeqpf , ZGEQPF )
#define F77_zgeqr2 AOCL_LAPACK_FUNC( zgeqr2 , ZGEQR2 )
#define F77_zgeqr2p AOCL_LAPACK_FUNC( zgeqr2p , ZGEQR2P )
#define F77_zgeqr AOCL_LAPACK_FUNC( zgeqr , ZGEQR )
#define F77_zgeqrf AOCL_LAPACK_FUNC( zgeqrf , ZGEQRF )
#define F77_zgeqrfp AOCL_LAPACK_FUNC( zgeqrfp , ZGEQRFP )
#define F77_zgeqrt2 AOCL_LAPACK_FUNC( zgeqrt2 , ZGEQRT2 )
#define F77_zgeqrt3 AOCL_LAPACK_FUNC( zgeqrt3 , ZGEQRT3 )
#define F77_zgeqrt AOCL_LAPACK_FUNC( zgeqrt , ZGEQRT )
#define F77_zgerfs AOCL_LAPACK_FUNC( zgerfs , ZGERFS )
#define F77_zgerfsx F77_FUNC( zgerfsx , ZGERFSX )
#define F77_zgerq2 AOCL_LAPACK_FUNC( zgerq2 , ZGERQ2 )
#define F77_zgerqf AOCL_LAPACK_FUNC( zgerqf , ZGERQF )
#define F77_zgesc2 AOCL_LAPACK_FUNC( zgesc2 , ZGESC2 )
#define F77_zgesdd AOCL_LAPACK_FUNC( zgesdd , ZGESDD )
#define F77_zgesvd AOCL_LAPACK_FUNC( zgesvd , ZGESVD )
#define F77_zgesvdq AOCL_LAPACK_FUNC( zgesvdq , ZGESVDQ )
#define F77_zgesvdx AOCL_LAPACK_FUNC( zgesvdx , ZGESVDX )
#define F77_zgesv AOCL_LAPACK_FUNC( zgesv , ZGESV )
#define F77_zgesvj AOCL_LAPACK_FUNC( zgesvj , ZGESVJ )
#define F77_zgesvx AOCL_LAPACK_FUNC( zgesvx , ZGESVX )
#define F77_zgesvxx F77_FUNC( zgesvxx , ZGESVXX )
#define F77_zgetc2 AOCL_LAPACK_FUNC( zgetc2 , ZGETC2 )
#define F77_zgetf2 AOCL_LAPACK_FUNC( zgetf2 , ZGETF2 )
#define F77_zgetrf2 AOCL_LAPACK_FUNC( zgetrf2 , ZGETRF2 )
#define F77_zgetrf AOCL_LAPACK_FUNC( zgetrf , ZGETRF )
#define F77_zgetri AOCL_LAPACK_FUNC( zgetri , ZGETRI )
#define F77_zgetrs AOCL_LAPACK_FUNC( zgetrs , ZGETRS )
#define F77_zgetsls AOCL_LAPACK_FUNC( zgetsls , ZGETSLS )
#define F77_zggbak AOCL_LAPACK_FUNC( zggbak , ZGGBAK )
#define F77_zggbal AOCL_LAPACK_FUNC( zggbal , ZGGBAL )
#define F77_zgges3 AOCL_LAPACK_FUNC( zgges3 , ZGGES3 )
#define F77_zgges AOCL_LAPACK_FUNC( zgges , ZGGES )
#define F77_zggesx AOCL_LAPACK_FUNC( zggesx , ZGGESX )
#define F77_zggev3 AOCL_LAPACK_FUNC( zggev3 , ZGGEV3 )
#define F77_zggev AOCL_LAPACK_FUNC( zggev , ZGGEV )
#define F77_zggevx AOCL_LAPACK_FUNC( zggevx , ZGGEVX )
#define F77_zggglm AOCL_LAPACK_FUNC( zggglm , ZGGGLM )
#define F77_zgghd3 AOCL_LAPACK_FUNC( zgghd3 , ZGGHD3 )
#define F77_zgghrd AOCL_LAPACK_FUNC( zgghrd , ZGGHRD )
#define F77_zgglse AOCL_LAPACK_FUNC( zgglse , ZGGLSE )
#define F77_zggqrf AOCL_LAPACK_FUNC( zggqrf , ZGGQRF )
#define F77_zggrqf AOCL_LAPACK_FUNC( zggrqf , ZGGRQF )
#define F77_zggsvd3 AOCL_LAPACK_FUNC( zggsvd3 , ZGGSVD3 )
#define F77_zggsvd AOCL_LAPACK_FUNC( zggsvd , ZGGSVD )
#define F77_zggsvp3 AOCL_LAPACK_FUNC( zggsvp3 , ZGGSVP3 )
#define F77_zggsvp AOCL_LAPACK_FUNC( zggsvp , ZGGSVP )
#define F77_zgsvj0 AOCL_LAPACK_FUNC( zgsvj0 , ZGSVJ0 )
#define F77_zgsvj1 AOCL_LAPACK_FUNC( zgsvj1 , ZGSVJ1 )
#define F77_zgtcon AOCL_LAPACK_FUNC( zgtcon , ZGTCON )
#define F77_zgtrfs AOCL_LAPACK_FUNC( zgtrfs , ZGTRFS )
#define F77_zgtsv AOCL_LAPACK_FUNC( zgtsv , ZGTSV )
#define F77_zgtsvx AOCL_LAPACK_FUNC( zgtsvx , ZGTSVX )
#define F77_zgttrf AOCL_LAPACK_FUNC( zgttrf , ZGTTRF )
#define F77_zgttrs AOCL_LAPACK_FUNC( zgttrs , ZGTTRS )
#define F77_zgtts2 AOCL_LAPACK_FUNC( zgtts2 , ZGTTS2 )
#define F77_zhb2st_kernels F77_FUNC( zhb2st_kernels , ZHB2ST_KERNELS )
#define F77_zhbev_2stage AOCL_LAPACK_FUNC( zhbev_2stage , ZHBEV_2STAGE )
#define F77_zhbevd_2stage AOCL_LAPACK_FUNC( zhbevd_2stage , ZHBEVD_2STAGE )
#define F77_zhbevd AOCL_LAPACK_FUNC( zhbevd , ZHBEVD )
#define F77_zhbev AOCL_LAPACK_FUNC( zhbev , ZHBEV )
#define F77_zhbevx_2stage AOCL_LAPACK_FUNC( zhbevx_2stage , ZHBEVX_2STAGE )
#define F77_zhbevx AOCL_LAPACK_FUNC( zhbevx , ZHBEVX )
#define F77_zhbgst AOCL_LAPACK_FUNC( zhbgst , ZHBGST )
#define F77_zhbgvd AOCL_LAPACK_FUNC( zhbgvd , ZHBGVD )
#define F77_zhbgv AOCL_LAPACK_FUNC( zhbgv , ZHBGV )
#define F77_zhbgvx AOCL_LAPACK_FUNC( zhbgvx , ZHBGVX )
#define F77_zhbtrd AOCL_LAPACK_FUNC( zhbtrd , ZHBTRD )
#define F77_zhecon_3 AOCL_LAPACK_FUNC( zhecon_3 , ZHECON_3 )
#define F77_zhecon AOCL_LAPACK_FUNC( zhecon , ZHECON )
#define F77_zhecon_rook AOCL_LAPACK_FUNC( zhecon_rook , ZHECON_ROOK )
#define F77_zheequb AOCL_LAPACK_FUNC( zheequb , ZHEEQUB )
#define F77_zheev_2stage AOCL_LAPACK_FUNC( zheev_2stage , ZHEEV_2STAGE )
#define F77_zheevd_2stage AOCL_LAPACK_FUNC( zheevd_2stage , ZHEEVD_2STAGE )
#define F77_zheevd AOCL_LAPACK_FUNC( zheevd , ZHEEVD )
#define F77_zheev  AOCL_LAPACK_FUNC( zheev  , ZHEEV  )
#define F77_zheev AOCL_LAPACK_FUNC( zheev , ZHEEV )
#define F77_zheevr_2stage AOCL_LAPACK_FUNC( zheevr_2stage , ZHEEVR_2STAGE )
#define F77_zheevr AOCL_LAPACK_FUNC( zheevr , ZHEEVR )
#define F77_zheevx_2stage AOCL_LAPACK_FUNC( zheevx_2stage , ZHEEVX_2STAGE )
#define F77_zheevx AOCL_LAPACK_FUNC( zheevx , ZHEEVX )
#define F77_zhegs2 AOCL_LAPACK_FUNC( zhegs2 , ZHEGS2 )
#define F77_zhegst AOCL_LAPACK_FUNC( zhegst , ZHEGST )
#define F77_zhegv_2stage AOCL_LAPACK_FUNC( zhegv_2stage , ZHEGV_2STAGE )
#define F77_zhegvd AOCL_LAPACK_FUNC( zhegvd , ZHEGVD )
#define F77_zhegv AOCL_LAPACK_FUNC( zhegv , ZHEGV )
#define F77_zhegvx AOCL_LAPACK_FUNC( zhegvx , ZHEGVX )
#define F77_zherfs AOCL_LAPACK_FUNC( zherfs , ZHERFS )
#define F77_zherfsx F77_FUNC( zherfsx , ZHERFSX )
#define F77_zhesv_aa_2stage AOCL_LAPACK_FUNC( zhesv_aa_2stage , ZHESV_AA_2STAGE )
#define F77_zhesv_aa AOCL_LAPACK_FUNC( zhesv_aa , ZHESV_AA )
#define F77_zhesv AOCL_LAPACK_FUNC( zhesv , ZHESV )
#define F77_zhesv_rk AOCL_LAPACK_FUNC( zhesv_rk , ZHESV_RK )
#define F77_zhesv_rook AOCL_LAPACK_FUNC( zhesv_rook , ZHESV_ROOK )
#define F77_zhesvx AOCL_LAPACK_FUNC( zhesvx , ZHESVX )
#define F77_zhesvxx F77_FUNC( zhesvxx , ZHESVXX )
#define F77_zheswapr AOCL_LAPACK_FUNC( zheswapr , ZHESWAPR )
#define F77_zhetd2 AOCL_LAPACK_FUNC( zhetd2 , ZHETD2 )
#define F77_zhetf2 AOCL_LAPACK_FUNC( zhetf2 , ZHETF2 )
#define F77_zhetf2_rk AOCL_LAPACK_FUNC( zhetf2_rk , ZHETF2_RK )
#define F77_zhetf2_rook AOCL_LAPACK_FUNC( zhetf2_rook , ZHETF2_ROOK )
#define F77_zhetrd_2stage AOCL_LAPACK_FUNC( zhetrd_2stage , ZHETRD_2STAGE )
#define F77_zhetrd AOCL_LAPACK_FUNC( zhetrd , ZHETRD )
#define F77_zhetrd_hb2st AOCL_LAPACK_FUNC( zhetrd_hb2st , ZHETRD_HB2ST )
#define F77_zhetrd_he2hb AOCL_LAPACK_FUNC( zhetrd_he2hb , ZHETRD_HE2HB )
#define F77_zhetrf_aa_2stage AOCL_LAPACK_FUNC( zhetrf_aa_2stage , ZHETRF_AA_2STAGE )
#define F77_zhetrf_aa AOCL_LAPACK_FUNC( zhetrf_aa , ZHETRF_AA )
#define F77_zhetrf AOCL_LAPACK_FUNC( zhetrf , ZHETRF )
#define F77_zhetrf_rk AOCL_LAPACK_FUNC( zhetrf_rk , ZHETRF_RK )
#define F77_zhetrf_rook AOCL_LAPACK_FUNC( zhetrf_rook , ZHETRF_ROOK )
#define F77_zhetri2 AOCL_LAPACK_FUNC( zhetri2 , ZHETRI2 )
#define F77_zhetri2x AOCL_LAPACK_FUNC( zhetri2x , ZHETRI2X )
#define F77_zhetri_3 AOCL_LAPACK_FUNC( zhetri_3 , ZHETRI_3 )
#define F77_zhetri_3x AOCL_LAPACK_FUNC( zhetri_3x , ZHETRI_3X )
#define F77_zhetri AOCL_LAPACK_FUNC( zhetri , ZHETRI )
#define F77_zhetri_rook AOCL_LAPACK_FUNC( zhetri_rook , ZHETRI_ROOK )
#define F77_zhetrs2 AOCL_LAPACK_FUNC( zhetrs2 , ZHETRS2 )
#define F77_zhetrs_3 AOCL_LAPACK_FUNC( zhetrs_3 , ZHETRS_3 )
#define F77_zhetrs_aa_2stage AOCL_LAPACK_FUNC( zhetrs_aa_2stage , ZHETRS_AA_2STAGE )
#define F77_zhetrs_aa AOCL_LAPACK_FUNC( zhetrs_aa , ZHETRS_AA )
#define F77_zhetrs AOCL_LAPACK_FUNC( zhetrs , ZHETRS )
#define F77_zhetrs_rook AOCL_LAPACK_FUNC( zhetrs_rook , ZHETRS_ROOK )
#define F77_zhfrk AOCL_LAPACK_FUNC( zhfrk , ZHFRK )
#define F77_zhgeqz AOCL_LAPACK_FUNC( zhgeqz , ZHGEQZ )
#define F77_zhpcon AOCL_LAPACK_FUNC( zhpcon , ZHPCON )
#define F77_zhpevd AOCL_LAPACK_FUNC( zhpevd , ZHPEVD )
#define F77_zhpev AOCL_LAPACK_FUNC( zhpev , ZHPEV )
#define F77_zhpevx AOCL_LAPACK_FUNC( zhpevx , ZHPEVX )
#define F77_zhpgst AOCL_LAPACK_FUNC( zhpgst , ZHPGST )
#define F77_zhpgvd AOCL_LAPACK_FUNC( zhpgvd , ZHPGVD )
#define F77_zhpgv AOCL_LAPACK_FUNC( zhpgv , ZHPGV )
#define F77_zhpgvx AOCL_LAPACK_FUNC( zhpgvx , ZHPGVX )
#define F77_zhprfs AOCL_LAPACK_FUNC( zhprfs , ZHPRFS )
#define F77_zhpsv AOCL_LAPACK_FUNC( zhpsv , ZHPSV )
#define F77_zhpsvx AOCL_LAPACK_FUNC( zhpsvx , ZHPSVX )
#define F77_zhptrd AOCL_LAPACK_FUNC( zhptrd , ZHPTRD )
#define F77_zhptrf AOCL_LAPACK_FUNC( zhptrf , ZHPTRF )
#define F77_zhptri AOCL_LAPACK_FUNC( zhptri , ZHPTRI )
#define F77_zhptrs AOCL_LAPACK_FUNC( zhptrs , ZHPTRS )
#define F77_zhsein AOCL_LAPACK_FUNC( zhsein , ZHSEIN )
#define F77_zhseqr AOCL_LAPACK_FUNC( zhseqr , ZHSEQR )
#define F77_zlabrd AOCL_LAPACK_FUNC( zlabrd , ZLABRD )
#define F77_zlacgv AOCL_LAPACK_FUNC( zlacgv , ZLACGV )
#define F77_zlacn2 AOCL_LAPACK_FUNC( zlacn2 , ZLACN2 )
#define F77_zlacon AOCL_LAPACK_FUNC( zlacon , ZLACON )
#define F77_zlacp2 AOCL_LAPACK_FUNC( zlacp2 , ZLACP2 )
#define F77_zlacpy AOCL_LAPACK_FUNC( zlacpy , ZLACPY )
#define F77_zlacrm AOCL_LAPACK_FUNC( zlacrm , ZLACRM )
#define F77_zlacrt AOCL_LAPACK_FUNC( zlacrt , ZLACRT )
#define F77_zladiv F77_FUNC( zladiv , ZLADIV )
#define F77_zlaed0 AOCL_LAPACK_FUNC( zlaed0 , ZLAED0 )
#define F77_zlaed7 AOCL_LAPACK_FUNC( zlaed7 , ZLAED7 )
#define F77_zlaed8 AOCL_LAPACK_FUNC( zlaed8 , ZLAED8 )
#define F77_zlaein AOCL_LAPACK_FUNC( zlaein , ZLAEIN )
#define F77_zlaesy F77_FUNC( zlaesy , ZLAESY )
#define F77_zlaev2 F77_FUNC( zlaev2 , ZLAEV2 )
#define F77_zlag2c AOCL_LAPACK_FUNC( zlag2c , ZLAG2C )
#define F77_zla_gbamv F77_FUNC( zla_gbamv , ZLA_GBAMV )
#define F77_zla_gbrcond_c F77_FUNC( zla_gbrcond_c , ZLA_GBRCOND_C )
#define F77_zla_gbrcond_x F77_FUNC( zla_gbrcond_x , ZLA_GBRCOND_X )
#define F77_zla_gbrfsx_extended F77_FUNC( zla_gbrfsx_extended , ZLA_GBRFSX_EXTENDED )
#define F77_zla_gbrpvgrw F77_FUNC( zla_gbrpvgrw , ZLA_GBRPVGRW )
#define F77_zla_geamv F77_FUNC( zla_geamv , ZLA_GEAMV )
#define F77_zla_gercond_c F77_FUNC( zla_gercond_c , ZLA_GERCOND_C )
#define F77_zla_gercond_x F77_FUNC( zla_gercond_x , ZLA_GERCOND_X )
#define F77_zla_gerfsx_extended F77_FUNC( zla_gerfsx_extended , ZLA_GERFSX_EXTENDED )
#define F77_zla_gerpvgrw F77_FUNC( zla_gerpvgrw , ZLA_GERPVGRW )
#define F77_zlags2 F77_FUNC( zlags2 , ZLAGS2 )
#define F77_zlagtm AOCL_LAPACK_FUNC( zlagtm , ZLAGTM )
#define F77_zla_heamv F77_FUNC( zla_heamv , ZLA_HEAMV )
#define F77_zlahef_aa AOCL_LAPACK_FUNC( zlahef_aa , ZLAHEF_AA )
#define F77_zlahef AOCL_LAPACK_FUNC( zlahef , ZLAHEF )
#define F77_zlahef_rk AOCL_LAPACK_FUNC( zlahef_rk , ZLAHEF_RK )
#define F77_zlahef_rook AOCL_LAPACK_FUNC( zlahef_rook , ZLAHEF_ROOK )
#define F77_zla_hercond_c F77_FUNC( zla_hercond_c , ZLA_HERCOND_C )
#define F77_zla_hercond_x F77_FUNC( zla_hercond_x , ZLA_HERCOND_X )
#define F77_zla_herfsx_extended F77_FUNC( zla_herfsx_extended , ZLA_HERFSX_EXTENDED )
#define F77_zla_herpvgrw F77_FUNC( zla_herpvgrw , ZLA_HERPVGRW )
#define F77_zlahqr AOCL_LAPACK_FUNC( zlahqr , ZLAHQR )
#define F77_zlahr2 AOCL_LAPACK_FUNC( zlahr2 , ZLAHR2 )
#define F77_zlahrd AOCL_LAPACK_FUNC( zlahrd , ZLAHRD )
#define F77_zlaic1 AOCL_LAPACK_FUNC( zlaic1 , ZLAIC1 )
#define F77_zla_lin_berr F77_FUNC( zla_lin_berr , ZLA_LIN_BERR )
#define F77_zlals0 AOCL_LAPACK_FUNC( zlals0 , ZLALS0 )
#define F77_zlalsa AOCL_LAPACK_FUNC( zlalsa , ZLALSA )
#define F77_zlalsd AOCL_LAPACK_FUNC( zlalsd , ZLALSD )
#define F77_zlamswlq AOCL_LAPACK_FUNC( zlamswlq , ZLAMSWLQ )
#define F77_zlamtsqr AOCL_LAPACK_FUNC( zlamtsqr , ZLAMTSQR )
#define F77_zlangb AOCL_LAPACK_FUNC( zlangb , ZLANGB )
#define F77_zlange AOCL_LAPACK_FUNC( zlange , ZLANGE )
#define F77_zlangt AOCL_LAPACK_FUNC( zlangt , ZLANGT )
#define F77_zlanhb AOCL_LAPACK_FUNC( zlanhb , ZLANHB )
#define F77_zlanhe AOCL_LAPACK_FUNC( zlanhe , ZLANHE )
#define F77_zlanhf AOCL_LAPACK_FUNC( zlanhf , ZLANHF )
#define F77_zlanhp AOCL_LAPACK_FUNC( zlanhp , ZLANHP )
#define F77_zlanhs AOCL_LAPACK_FUNC( zlanhs , ZLANHS )
#define F77_zlanht AOCL_LAPACK_FUNC( zlanht , ZLANHT )
#define F77_zlansb AOCL_LAPACK_FUNC( zlansb , ZLANSB )
#define F77_zlansp AOCL_LAPACK_FUNC( zlansp , ZLANSP )
#define F77_zlansy AOCL_LAPACK_FUNC( zlansy , ZLANSY )
#define F77_zlantb AOCL_LAPACK_FUNC( zlantb , ZLANTB )
#define F77_zlantp AOCL_LAPACK_FUNC( zlantp , ZLANTP )
#define F77_zlantr AOCL_LAPACK_FUNC( zlantr , ZLANTR )
#define F77_zlapll AOCL_LAPACK_FUNC( zlapll , ZLAPLL )
#define F77_zlapmr AOCL_LAPACK_FUNC( zlapmr , ZLAPMR )
#define F77_zlapmt AOCL_LAPACK_FUNC( zlapmt , ZLAPMT )
#define F77_zla_porcond_c F77_FUNC( zla_porcond_c , ZLA_PORCOND_C )
#define F77_zla_porcond_x F77_FUNC( zla_porcond_x , ZLA_PORCOND_X )
#define F77_zla_porfsx_extended F77_FUNC( zla_porfsx_extended , ZLA_PORFSX_EXTENDED )
#define F77_zla_porpvgrw F77_FUNC( zla_porpvgrw , ZLA_PORPVGRW )
#define F77_zlaqgb AOCL_LAPACK_FUNC( zlaqgb , ZLAQGB )
#define F77_zlaqge AOCL_LAPACK_FUNC( zlaqge , ZLAQGE )
#define F77_zlaqhb AOCL_LAPACK_FUNC( zlaqhb , ZLAQHB )
#define F77_zlaqhe AOCL_LAPACK_FUNC( zlaqhe , ZLAQHE )
#define F77_zlaqhp AOCL_LAPACK_FUNC( zlaqhp , ZLAQHP )
#define F77_zlaqp2 AOCL_LAPACK_FUNC( zlaqp2 , ZLAQP2 )
#define F77_zlaqps AOCL_LAPACK_FUNC( zlaqps , ZLAQPS )
#define F77_zlaqr0 AOCL_LAPACK_FUNC( zlaqr0 , ZLAQR0 )
#define F77_zlaqr1 AOCL_LAPACK_FUNC( zlaqr1 , ZLAQR1 )
#define F77_zlaqr2 AOCL_LAPACK_FUNC( zlaqr2 , ZLAQR2 )
#define F77_zlaqr3 AOCL_LAPACK_FUNC( zlaqr3 , ZLAQR3 )
#define F77_zlaqr4 AOCL_LAPACK_FUNC( zlaqr4 , ZLAQR4 )
#define F77_zlaqr5 AOCL_LAPACK_FUNC( zlaqr5 , ZLAQR5 )
#define F77_zlaqsb AOCL_LAPACK_FUNC( zlaqsb , ZLAQSB )
#define F77_zlaqsp AOCL_LAPACK_FUNC( zlaqsp , ZLAQSP )
#define F77_zlaqsy AOCL_LAPACK_FUNC( zlaqsy , ZLAQSY )
#define F77_zlar1v AOCL_LAPACK_FUNC( zlar1v , ZLAR1V )
#define F77_zlar2v AOCL_LAPACK_FUNC( zlar2v , ZLAR2V )
#define F77_zlarcm AOCL_LAPACK_FUNC( zlarcm , ZLARCM )
#define F77_zlarfb AOCL_LAPACK_FUNC( zlarfb , ZLARFB )
#define F77_zlarf AOCL_LAPACK_FUNC( zlarf , ZLARF )
#define F77_zlarfg AOCL_LAPACK_FUNC( zlarfg , ZLARFG )
#define F77_zlarfgp AOCL_LAPACK_FUNC( zlarfgp , ZLARFGP )
#define F77_zlarft AOCL_LAPACK_FUNC( zlarft , ZLARFT )
#define F77_zlarfx AOCL_LAPACK_FUNC( zlarfx , ZLARFX )
#define F77_zlarfy AOCL_LAPACK_FUNC( zlarfy , ZLARFY )
#define F77_zlargv AOCL_LAPACK_FUNC( zlargv , ZLARGV )
#define F77_zlarnv AOCL_LAPACK_FUNC( zlarnv , ZLARNV )
#define F77_zlarrv AOCL_LAPACK_FUNC( zlarrv , ZLARRV )
#define F77_zlarscl2 AOCL_LAPACK_FUNC( zlarscl2 , ZLARSCL2 )
#define F77_zlartg F77_FUNC( zlartg , ZLARTG )
#define F77_zlartv AOCL_LAPACK_FUNC( zlartv , ZLARTV )
#define F77_zlarzb AOCL_LAPACK_FUNC( zlarzb , ZLARZB )
#define F77_zlarz AOCL_LAPACK_FUNC( zlarz , ZLARZ )
#define F77_zlarzt AOCL_LAPACK_FUNC( zlarzt , ZLARZT )
#define F77_zlascl2 AOCL_LAPACK_FUNC( zlascl2 , ZLASCL2 )
#define F77_zlascl AOCL_LAPACK_FUNC( zlascl , ZLASCL )
#define F77_zlaset AOCL_LAPACK_FUNC( zlaset , ZLASET )
#define F77_zlasr AOCL_LAPACK_FUNC( zlasr , ZLASR )
#define F77_zlassq AOCL_LAPACK_FUNC( zlassq , ZLASSQ )
#define F77_zlaswlq AOCL_LAPACK_FUNC( zlaswlq , ZLASWLQ )
#define F77_zlaswp AOCL_LAPACK_FUNC( zlaswp , ZLASWP )
#define F77_zla_syamv F77_FUNC( zla_syamv , ZLA_SYAMV )
#define F77_zlasyf_aa AOCL_LAPACK_FUNC( zlasyf_aa , ZLASYF_AA )
#define F77_zlasyf AOCL_LAPACK_FUNC( zlasyf , ZLASYF )
#define F77_zlasyf_rk AOCL_LAPACK_FUNC( zlasyf_rk , ZLASYF_RK )
#define F77_zlasyf_rook AOCL_LAPACK_FUNC( zlasyf_rook , ZLASYF_ROOK )
#define F77_zla_syrcond_c F77_FUNC( zla_syrcond_c , ZLA_SYRCOND_C )
#define F77_zla_syrcond_x F77_FUNC( zla_syrcond_x , ZLA_SYRCOND_X )
#define F77_zla_syrfsx_extended F77_FUNC( zla_syrfsx_extended , ZLA_SYRFSX_EXTENDED )
#define F77_zla_syrpvgrw F77_FUNC( zla_syrpvgrw , ZLA_SYRPVGRW )
#define F77_zlat2c AOCL_LAPACK_FUNC( zlat2c , ZLAT2C )
#define F77_zlatbs AOCL_LAPACK_FUNC( zlatbs , ZLATBS )
#define F77_zlatdf AOCL_LAPACK_FUNC( zlatdf , ZLATDF )
#define F77_zlatps AOCL_LAPACK_FUNC( zlatps , ZLATPS )
#define F77_zlatrd AOCL_LAPACK_FUNC( zlatrd , ZLATRD )
#define F77_zlatrs AOCL_LAPACK_FUNC( zlatrs , ZLATRS )
#define F77_zlatrz AOCL_LAPACK_FUNC( zlatrz , ZLATRZ )
#define F77_zlatsqr AOCL_LAPACK_FUNC( zlatsqr , ZLATSQR )
#define F77_zlatzm AOCL_LAPACK_FUNC( zlatzm , ZLATZM )
#define F77_zlaunhr_col_getrfnp2 F77_FUNC( zlaunhr_col_getrfnp2 , ZLAUNHR_COL_GETRFNP2 )
#define F77_zlaunhr_col_getrfnp F77_FUNC( zlaunhr_col_getrfnp , ZLAUNHR_COL_GETRFNP )
#define F77_zlauu2 AOCL_LAPACK_FUNC( zlauu2 , ZLAUU2 )
#define F77_zlauum AOCL_LAPACK_FUNC( zlauum , ZLAUUM )
#define F77_zla_wwaddw F77_FUNC( zla_wwaddw , ZLA_WWADDW )
#define F77_zpbcon AOCL_LAPACK_FUNC( zpbcon , ZPBCON )
#define F77_zpbequ AOCL_LAPACK_FUNC( zpbequ , ZPBEQU )
#define F77_zpbrfs AOCL_LAPACK_FUNC( zpbrfs , ZPBRFS )
#define F77_zpbstf AOCL_LAPACK_FUNC( zpbstf , ZPBSTF )
#define F77_zpbsv AOCL_LAPACK_FUNC( zpbsv , ZPBSV )
#define F77_zpbsvx AOCL_LAPACK_FUNC( zpbsvx , ZPBSVX )
#define F77_zpbtf2 AOCL_LAPACK_FUNC( zpbtf2 , ZPBTF2 )
#define F77_zpbtrf AOCL_LAPACK_FUNC( zpbtrf , ZPBTRF )
#define F77_zpbtrs AOCL_LAPACK_FUNC( zpbtrs , ZPBTRS )
#define F77_zpftrf AOCL_LAPACK_FUNC( zpftrf , ZPFTRF )
#define F77_zpftri AOCL_LAPACK_FUNC( zpftri , ZPFTRI )
#define F77_zpftrs AOCL_LAPACK_FUNC( zpftrs , ZPFTRS )
#define F77_zpocon AOCL_LAPACK_FUNC( zpocon , ZPOCON )
#define F77_zpoequb AOCL_LAPACK_FUNC( zpoequb , ZPOEQUB )
#define F77_zpoequ AOCL_LAPACK_FUNC( zpoequ , ZPOEQU )
#define F77_zporfs AOCL_LAPACK_FUNC( zporfs , ZPORFS )
#define F77_zporfsx F77_FUNC( zporfsx , ZPORFSX )
#define F77_zposv AOCL_LAPACK_FUNC( zposv , ZPOSV )
#define F77_zposvx AOCL_LAPACK_FUNC( zposvx , ZPOSVX )
#define F77_zposvxx F77_FUNC( zposvxx , ZPOSVXX )
#define F77_zpotf2 AOCL_LAPACK_FUNC( zpotf2 , ZPOTF2 )
#define F77_zpotrf2 AOCL_LAPACK_FUNC( zpotrf2 , ZPOTRF2 )
#define F77_zpotrf AOCL_LAPACK_FUNC( zpotrf , ZPOTRF )
#define F77_zpotri AOCL_LAPACK_FUNC( zpotri , ZPOTRI )
#define F77_zpotrs AOCL_LAPACK_FUNC( zpotrs , ZPOTRS )
#define F77_zppcon AOCL_LAPACK_FUNC( zppcon , ZPPCON )
#define F77_zppequ AOCL_LAPACK_FUNC( zppequ , ZPPEQU )
#define F77_zpprfs AOCL_LAPACK_FUNC( zpprfs , ZPPRFS )
#define F77_zppsv AOCL_LAPACK_FUNC( zppsv , ZPPSV )
#define F77_zppsvx AOCL_LAPACK_FUNC( zppsvx , ZPPSVX )
#define F77_zpptrf AOCL_LAPACK_FUNC( zpptrf , ZPPTRF )
#define F77_zpptri AOCL_LAPACK_FUNC( zpptri , ZPPTRI )
#define F77_zpptrs AOCL_LAPACK_FUNC( zpptrs , ZPPTRS )
#define F77_zpstf2 AOCL_LAPACK_FUNC( zpstf2 , ZPSTF2 )
#define F77_zpstrf AOCL_LAPACK_FUNC( zpstrf , ZPSTRF )
#define F77_zptcon AOCL_LAPACK_FUNC( zptcon , ZPTCON )
#define F77_zpteqr AOCL_LAPACK_FUNC( zpteqr , ZPTEQR )
#define F77_zptrfs AOCL_LAPACK_FUNC( zptrfs , ZPTRFS )
#define F77_zptsv AOCL_LAPACK_FUNC( zptsv , ZPTSV )
#define F77_zptsvx AOCL_LAPACK_FUNC( zptsvx , ZPTSVX )
#define F77_zpttrf AOCL_LAPACK_FUNC( zpttrf , ZPTTRF )
#define F77_zpttrs AOCL_LAPACK_FUNC( zpttrs , ZPTTRS )
#define F77_zptts2 AOCL_LAPACK_FUNC( zptts2 , ZPTTS2 )
#define F77_zrot AOCL_LAPACK_FUNC( zrot , ZROT )
#define F77_zspcon AOCL_LAPACK_FUNC( zspcon , ZSPCON )
#define F77_zspmv AOCL_LAPACK_FUNC( zspmv , ZSPMV )
#define F77_zspr AOCL_LAPACK_FUNC( zspr , ZSPR )
#define F77_zsprfs AOCL_LAPACK_FUNC( zsprfs , ZSPRFS )
#define F77_zspsv AOCL_LAPACK_FUNC( zspsv , ZSPSV )
#define F77_zspsvx AOCL_LAPACK_FUNC( zspsvx , ZSPSVX )
#define F77_zsptrf AOCL_LAPACK_FUNC( zsptrf , ZSPTRF )
#define F77_zsptri AOCL_LAPACK_FUNC( zsptri , ZSPTRI )
#define F77_zsptrs AOCL_LAPACK_FUNC( zsptrs , ZSPTRS )
#define F77_zstedc AOCL_LAPACK_FUNC( zstedc , ZSTEDC )
#define F77_zstegr AOCL_LAPACK_FUNC( zstegr , ZSTEGR )
#define F77_zstein AOCL_LAPACK_FUNC( zstein , ZSTEIN )
#define F77_zstemr AOCL_LAPACK_FUNC( zstemr , ZSTEMR )
#define F77_zsteqr AOCL_LAPACK_FUNC( zsteqr , ZSTEQR )
#define F77_zsycon_3 AOCL_LAPACK_FUNC( zsycon_3 , ZSYCON_3 )
#define F77_zsycon AOCL_LAPACK_FUNC( zsycon , ZSYCON )
#define F77_zsycon_rook AOCL_LAPACK_FUNC( zsycon_rook , ZSYCON_ROOK )
#define F77_zsyconv AOCL_LAPACK_FUNC( zsyconv , ZSYCONV )
#define F77_zsyconvf AOCL_LAPACK_FUNC( zsyconvf , ZSYCONVF )
#define F77_zsyconvf_rook AOCL_LAPACK_FUNC( zsyconvf_rook , ZSYCONVF_ROOK )
#define F77_zsyequb AOCL_LAPACK_FUNC( zsyequb , ZSYEQUB )
#define F77_zsymv AOCL_LAPACK_FUNC( zsymv , ZSYMV )
#define F77_zsyr AOCL_LAPACK_FUNC( zsyr , ZSYR )
#define F77_zsyrfs AOCL_LAPACK_FUNC( zsyrfs , ZSYRFS )
#define F77_zsyrfsx F77_FUNC( zsyrfsx , ZSYRFSX )
#define F77_zsysv_aa_2stage AOCL_LAPACK_FUNC( zsysv_aa_2stage , ZSYSV_AA_2STAGE )
#define F77_zsysv_aa AOCL_LAPACK_FUNC( zsysv_aa , ZSYSV_AA )
#define F77_zsysv AOCL_LAPACK_FUNC( zsysv , ZSYSV )
#define F77_zsysv_rk AOCL_LAPACK_FUNC( zsysv_rk , ZSYSV_RK )
#define F77_zsysv_rook AOCL_LAPACK_FUNC( zsysv_rook , ZSYSV_ROOK )
#define F77_zsysvx AOCL_LAPACK_FUNC( zsysvx , ZSYSVX )
#define F77_zsysvxx F77_FUNC( zsysvxx , ZSYSVXX )
#define F77_zsyswapr AOCL_LAPACK_FUNC( zsyswapr , ZSYSWAPR )
#define F77_zsytf2 AOCL_LAPACK_FUNC( zsytf2 , ZSYTF2 )
#define F77_zsytf2_rk AOCL_LAPACK_FUNC( zsytf2_rk , ZSYTF2_RK )
#define F77_zsytf2_rook AOCL_LAPACK_FUNC( zsytf2_rook , ZSYTF2_ROOK )
#define F77_zsytrf_aa_2stage AOCL_LAPACK_FUNC( zsytrf_aa_2stage , ZSYTRF_AA_2STAGE )
#define F77_zsytrf_aa AOCL_LAPACK_FUNC( zsytrf_aa , ZSYTRF_AA )
#define F77_zsytrf AOCL_LAPACK_FUNC( zsytrf , ZSYTRF )
#define F77_zsytrf_rk AOCL_LAPACK_FUNC( zsytrf_rk , ZSYTRF_RK )
#define F77_zsytrf_rook AOCL_LAPACK_FUNC( zsytrf_rook , ZSYTRF_ROOK )
#define F77_zsytri2 AOCL_LAPACK_FUNC( zsytri2 , ZSYTRI2 )
#define F77_zsytri2x AOCL_LAPACK_FUNC( zsytri2x , ZSYTRI2X )
#define F77_zsytri_3 AOCL_LAPACK_FUNC( zsytri_3 , ZSYTRI_3 )
#define F77_zsytri_3x AOCL_LAPACK_FUNC( zsytri_3x , ZSYTRI_3X )
#define F77_zsytri AOCL_LAPACK_FUNC( zsytri , ZSYTRI )
#define F77_zsytri_rook AOCL_LAPACK_FUNC( zsytri_rook , ZSYTRI_ROOK )
#define F77_zsytrs2 AOCL_LAPACK_FUNC( zsytrs2 , ZSYTRS2 )
#define F77_zsytrs_3 AOCL_LAPACK_FUNC( zsytrs_3 , ZSYTRS_3 )
#define F77_zsytrs_aa_2stage AOCL_LAPACK_FUNC( zsytrs_aa_2stage , ZSYTRS_AA_2STAGE )
#define F77_zsytrs_aa AOCL_LAPACK_FUNC( zsytrs_aa , ZSYTRS_AA )
#define F77_zsytrs AOCL_LAPACK_FUNC( zsytrs , ZSYTRS )
#define F77_zsytrs_rook AOCL_LAPACK_FUNC( zsytrs_rook , ZSYTRS_ROOK )
#define F77_ztbcon AOCL_LAPACK_FUNC( ztbcon , ZTBCON )
#define F77_ztbrfs AOCL_LAPACK_FUNC( ztbrfs , ZTBRFS )
#define F77_ztbtrs AOCL_LAPACK_FUNC( ztbtrs , ZTBTRS )
#define F77_ztfsm AOCL_LAPACK_FUNC( ztfsm , ZTFSM )
#define F77_ztftri AOCL_LAPACK_FUNC( ztftri , ZTFTRI )
#define F77_ztfttp AOCL_LAPACK_FUNC( ztfttp , ZTFTTP )
#define F77_ztfttr AOCL_LAPACK_FUNC( ztfttr , ZTFTTR )
#define F77_ztgevc AOCL_LAPACK_FUNC( ztgevc , ZTGEVC )
#define F77_ztgex2 AOCL_LAPACK_FUNC( ztgex2 , ZTGEX2 )
#define F77_ztgexc AOCL_LAPACK_FUNC( ztgexc , ZTGEXC )
#define F77_ztgsen AOCL_LAPACK_FUNC( ztgsen , ZTGSEN )
#define F77_ztgsja AOCL_LAPACK_FUNC( ztgsja , ZTGSJA )
#define F77_ztgsna AOCL_LAPACK_FUNC( ztgsna , ZTGSNA )
#define F77_ztgsy2 AOCL_LAPACK_FUNC( ztgsy2 , ZTGSY2 )
#define F77_ztgsyl AOCL_LAPACK_FUNC( ztgsyl , ZTGSYL )
#define F77_ztpcon AOCL_LAPACK_FUNC( ztpcon , ZTPCON )
#define F77_ztplqt2 AOCL_LAPACK_FUNC( ztplqt2 , ZTPLQT2 )
#define F77_ztplqt AOCL_LAPACK_FUNC( ztplqt , ZTPLQT )
#define F77_ztpmlqt AOCL_LAPACK_FUNC( ztpmlqt , ZTPMLQT )
#define F77_ztpmqrt AOCL_LAPACK_FUNC( ztpmqrt , ZTPMQRT )
#define F77_ztpqrt2 AOCL_LAPACK_FUNC( ztpqrt2 , ZTPQRT2 )
#define F77_ztpqrt AOCL_LAPACK_FUNC( ztpqrt , ZTPQRT )
#define F77_ztprfb AOCL_LAPACK_FUNC( ztprfb , ZTPRFB )
#define F77_ztprfs AOCL_LAPACK_FUNC( ztprfs , ZTPRFS )
#define F77_ztptri AOCL_LAPACK_FUNC( ztptri , ZTPTRI )
#define F77_ztptrs AOCL_LAPACK_FUNC( ztptrs , ZTPTRS )
#define F77_ztpttf AOCL_LAPACK_FUNC( ztpttf , ZTPTTF )
#define F77_ztpttr AOCL_LAPACK_FUNC( ztpttr , ZTPTTR )
#define F77_ztrcon AOCL_LAPACK_FUNC( ztrcon , ZTRCON )
#define F77_ztrevc3 AOCL_LAPACK_FUNC( ztrevc3 , ZTREVC3 )
#define F77_ztrevc AOCL_LAPACK_FUNC( ztrevc , ZTREVC )
#define F77_ztrexc AOCL_LAPACK_FUNC( ztrexc , ZTREXC )
#define F77_ztrrfs AOCL_LAPACK_FUNC( ztrrfs , ZTRRFS )
#define F77_ztrsen AOCL_LAPACK_FUNC( ztrsen , ZTRSEN )
#define F77_ztrsna AOCL_LAPACK_FUNC( ztrsna , ZTRSNA )
#define F77_ztrsyl AOCL_LAPACK_FUNC( ztrsyl , ZTRSYL )
#define F77_ztrti2 AOCL_LAPACK_FUNC( ztrti2 , ZTRTI2 )
#define F77_ztrtri AOCL_LAPACK_FUNC( ztrtri , ZTRTRI )
#define F77_ztrtrs AOCL_LAPACK_FUNC( ztrtrs , ZTRTRS )
#define F77_ztrttf AOCL_LAPACK_FUNC( ztrttf , ZTRTTF )
#define F77_ztrttp AOCL_LAPACK_FUNC( ztrttp , ZTRTTP )
#define F77_ztzrqf AOCL_LAPACK_FUNC( ztzrqf , ZTZRQF )
#define F77_ztzrzf AOCL_LAPACK_FUNC( ztzrzf , ZTZRZF )
#define F77_zunbdb1 AOCL_LAPACK_FUNC( zunbdb1 , ZUNBDB1 )
#define F77_zunbdb2 AOCL_LAPACK_FUNC( zunbdb2 , ZUNBDB2 )
#define F77_zunbdb3 AOCL_LAPACK_FUNC( zunbdb3 , ZUNBDB3 )
#define F77_zunbdb4 AOCL_LAPACK_FUNC( zunbdb4 , ZUNBDB4 )
#define F77_zunbdb5 AOCL_LAPACK_FUNC( zunbdb5 , ZUNBDB5 )
#define F77_zunbdb6 AOCL_LAPACK_FUNC( zunbdb6 , ZUNBDB6 )
#define F77_zunbdb AOCL_LAPACK_FUNC( zunbdb , ZUNBDB )
#define F77_zuncsd2by1 AOCL_LAPACK_FUNC( zuncsd2by1 , ZUNCSD2BY1 )
#define F77_zuncsd AOCL_LAPACK_FUNC( zuncsd , ZUNCSD )
#define F77_zung2l AOCL_LAPACK_FUNC( zung2l , ZUNG2L )
#define F77_zung2r AOCL_LAPACK_FUNC( zung2r , ZUNG2R )
#define F77_zungbr AOCL_LAPACK_FUNC( zungbr , ZUNGBR )
#define F77_zunghr AOCL_LAPACK_FUNC( zunghr , ZUNGHR )
#define F77_zungl2 AOCL_LAPACK_FUNC( zungl2 , ZUNGL2 )
#define F77_zunglq AOCL_LAPACK_FUNC( zunglq , ZUNGLQ )
#define F77_zungql AOCL_LAPACK_FUNC( zungql , ZUNGQL )
#define F77_zungqr AOCL_LAPACK_FUNC( zungqr , ZUNGQR )
#define F77_zungr2 AOCL_LAPACK_FUNC( zungr2 , ZUNGR2 )
#define F77_zungrq AOCL_LAPACK_FUNC( zungrq , ZUNGRQ )
#define F77_zungtr AOCL_LAPACK_FUNC( zungtr , ZUNGTR )
#define F77_zungtsqr AOCL_LAPACK_FUNC( zungtsqr , ZUNGTSQR )
#define F77_zunhr_col F77_FUNC( zunhr_col , ZUNHR_COL )
#define F77_zunm22 AOCL_LAPACK_FUNC( zunm22 , ZUNM22 )
#define F77_zunm2l AOCL_LAPACK_FUNC( zunm2l , ZUNM2L )
#define F77_zunm2r AOCL_LAPACK_FUNC( zunm2r , ZUNM2R )
#define F77_zunmbr AOCL_LAPACK_FUNC( zunmbr , ZUNMBR )
#define F77_zunmhr AOCL_LAPACK_FUNC( zunmhr , ZUNMHR )
#define F77_zunml2 AOCL_LAPACK_FUNC( zunml2 , ZUNML2 )
#define F77_zunmlq AOCL_LAPACK_FUNC( zunmlq , ZUNMLQ )
#define F77_zunmql AOCL_LAPACK_FUNC( zunmql , ZUNMQL )
#define F77_zunmqr AOCL_LAPACK_FUNC( zunmqr , ZUNMQR )
#define F77_zunmr2 AOCL_LAPACK_FUNC( zunmr2 , ZUNMR2 )
#define F77_zunmr3 AOCL_LAPACK_FUNC( zunmr3 , ZUNMR3 )
#define F77_zunmrq AOCL_LAPACK_FUNC( zunmrq , ZUNMRQ )
#define F77_zunmrz AOCL_LAPACK_FUNC( zunmrz , ZUNMRZ )
#define F77_zunmtr AOCL_LAPACK_FUNC( zunmtr , ZUNMTR )
#define F77_zupgtr AOCL_LAPACK_FUNC( zupgtr , ZUPGTR )
#define F77_zupmtr AOCL_LAPACK_FUNC( zupmtr , ZUPMTR )
#define F77_sgetrfnp AOCL_LAPAKC_FUNC( sgetrfnp , SGETRFNP )
#define F77_dgetrfnp AOCL_LAPAKC_FUNC( dgetrfnp , DGETRFNP )
#define F77_cgetrfnp AOCL_LAPAKC_FUNC( cgetrfnp , CGETRFNP )
#define F77_zgetrfnp AOCL_LAPAKC_FUNC( zgetrfnp , ZGETRFNP )
#define F77_sspffrt2 AOCL_LAPAKC_FUNC( sspffrt2 , SSPFFRT2 )
#define F77_dspffrt2 AOCL_LAPAKC_FUNC( dspffrt2 , DSPFFRT2 )
#define F77_cspffrt2 AOCL_LAPAKC_FUNC( cspffrt2 , CSPFFRT2 )
#define F77_zspffrt2 AOCL_LAPAKC_FUNC( zspffrt2 , ZSPFFRT2 )
#define F77_sspffrtx AOCL_LAPAKC_FUNC( sspffrtx , SSPFFRTX )
#define F77_dspffrtx AOCL_LAPAKC_FUNC( dspffrtx , DSPFFRTX )
#define F77_cspffrtx AOCL_LAPAKC_FUNC( cspffrtx , CSPFFRTX )
#define F77_zspffrtx AOCL_LAPAKC_FUNC( zspffrtx , ZSPFFRTX )
#define F77_cgetrfnpi AOCL_LAPAKC_FUNC( cgetrfnpi , CGETRFNPI )
#define F77_dgetrfnpi AOCL_LAPAKC_FUNC( dgetrfnpi , DGETRFNPI )
#define F77_sgetrfnpi AOCL_LAPAKC_FUNC( sgetrfnpi , SGETRFNPI )
#define F77_zgetrfnpi AOCL_LAPAKC_FUNC( zgetrfnpi , ZGETRFNPI )

#define F77_sopgtr AOCL_LAPACK_FUNC( sopgtr, SOPGTR )
#define F77_dopgtr AOCL_LAPACK_FUNC( dopgtr, DOPGTR )
#define F77_sorcsd AOCL_LAPACK_FUNC( sorcsd, SORCSD )
#define F77_dorcsd AOCL_LAPACK_FUNC( dorcsd, DORCSD )
#define F77_sorcsd2by1 AOCL_LAPACK_FUNC( sorcsd2by1, SORCSD2BY1 )
#define F77_dorcsd2by1 AOCL_LAPACK_FUNC( dorcsd2by1, DORCSD2BY1 )
#define F77_sorghr AOCL_LAPACK_FUNC( sorghr, SORGHR )
#define F77_dorghr AOCL_LAPACK_FUNC( dorghr, DORGHR )
#define F77_sormhr AOCL_LAPACK_FUNC( sormhr, SORMHR )
#define F77_dormhr AOCL_LAPACK_FUNC( dormhr, DORMHR )

#define F77_sgedmd AOCL_LAPACK_FUNC( sgedmd, SGEDMD )
#define F77_dgedmd AOCL_LAPACK_FUNC( dgedmd, DGEDMD )
#define F77_cgedmd AOCL_LAPACK_FUNC( cgedmd, CGEDMD )
#define F77_zgedmd AOCL_LAPACK_FUNC( zgedmd, ZGEDMD )

#define F77_sgedmdq AOCL_LAPACK_FUNC( sgedmdq, SGEDMDQ )
#define F77_dgedmdq AOCL_LAPACK_FUNC( dgedmdq, DGEDMDQ )
#define F77_cgedmdq AOCL_LAPACK_FUNC( cgedmdq, CGEDMDQ )
#define F77_zgedmdq AOCL_LAPACK_FUNC( zgedmdq, ZGEDMDQ )



#define LAPACK_EXPORT_cgelst F77_FUNC( cgelst , CGELST )
#define LAPACK_EXPORT_clatrs3 F77_FUNC( clatrs3 , CLATRS3 )
#define LAPACK_EXPORT_ctrsyl3 F77_FUNC( ctrsyl3 , CTRSYL3 )
#define LAPACK_EXPORT_dlarmm F77_FUNC( dlarmm , DLARMM )
#define LAPACK_EXPORT_dgelst F77_FUNC( dgelst , DGELST )
#define LAPACK_EXPORT_dlatrs3 F77_FUNC( dlatrs3 , DLATRS3 )
#define LAPACK_EXPORT_dtrsyl3 F77_FUNC( dtrsyl3 , DTRSYL3 )
#define LAPACK_EXPORT_slarmm F77_FUNC( slarmm , SLARMM )
#define LAPACK_EXPORT_sgelst F77_FUNC( sgelst , SGELST )
#define LAPACK_EXPORT_slatrs3 F77_FUNC( slatrs3 , SLATRS3 )
#define LAPACK_EXPORT_strsyl3 F77_FUNC( strsyl3 , STRSYL3 )
#define LAPACK_EXPORT_zgelst F77_FUNC( zgelst , ZGELST )
#define LAPACK_EXPORT_zlatrs3 F77_FUNC( zlatrs3 , ZLATRS3 )
#define LAPACK_EXPORT_ztrsyl3 F77_FUNC( ztrsyl3 , ZTRSYL3 )
#define LAPACK_EXPORT_dlamch F77_FUNC( dlamch , DLAMCH )
#define LAPACK_EXPORT_dlamc3 F77_FUNC( dlamc3 , DLAMC3 )
#define LAPACK_EXPORT_dladiv2 F77_FUNC( dladiv2 , DLADIV2 )
#define LAPACK_EXPORT_dladiv1 F77_FUNC(dladiv1 , DLADIV1 )
#define LAPACK_EXPORT_sladiv2 F77_FUNC( sladiv2 , SLADIV2 )
#define LAPACK_EXPORT_sladiv1 F77_FUNC( sladiv1 , SLADIV1 )
#define LAPACK_EXPORT_iparmq F77_FUNC( iparmq, IPARMQ )
#define LAPACK_EXPORT_ilaenv F77_FUNC( ilaenv , ILAENV )
#define LAPACK_EXPORT_ieeeck F77_FUNC( ieeeck , IEEECK )
#define LAPACK_EXPORT_lsamen F77_FUNC( lsamen , LSAMEN )
#define LAPACK_EXPORT_slamc3 F77_FUNC( slamc3 , SLAMC3 )
#define LAPACK_EXPORT_slamch F77_FUNC( slamch, SLAMCH )
#define LAPACK_EXPORT_cgetsqrhrt F77_FUNC( cgetsqrhrt , CGETSQRHRT )
#define LAPACK_EXPORT_claqz0 F77_FUNC( claqz0 , CLAQZ0 )
#define LAPACK_EXPORT_claqz1 F77_FUNC( claqz1 , CLAQZ1 )
#define LAPACK_EXPORT_claqz2 F77_FUNC( claqz2 , CLAQZ2 )
#define LAPACK_EXPORT_claqz3 F77_FUNC( claqz3 , CLAQZ3 )
#define LAPACK_EXPORT_clarfb_gett F77_FUNC( clarfb_gett , CLARFB_GETT )
#define LAPACK_EXPORT_cungtsqr_row F77_FUNC( cungtsqr_row , CUNGTSQR_ROW )
#define LAPACK_EXPORT_dgetsqrhrt F77_FUNC( dgetsqrhrt , DGETSQRHRT )
#define LAPACK_EXPORT_dlaqz0 F77_FUNC( dlaqz0 , DLAQZ0 )
#define LAPACK_EXPORT_dlaqz1 F77_FUNC( dlaqz1 , DLAQZ1 )
#define LAPACK_EXPORT_dlaqz2 F77_FUNC( dlaqz2 , DLAQZ2 )
#define LAPACK_EXPORT_dlaqz3 F77_FUNC( dlaqz3 , DLAQZ3 )
#define LAPACK_EXPORT_dlaqz4 F77_FUNC( dlaqz4 , DLAQZ4 )
#define LAPACK_EXPORT_dlarfb_gett F77_FUNC( dlarfb_gett , DLARFB_GETT )
#define LAPACK_EXPORT_dorgtsqr_row F77_FUNC( dorgtsqr_row , DORGTSQR_ROW )
#define LAPACK_EXPORT_sgetsqrhrt F77_FUNC( sgetsqrhrt , SGETSQRHRT )
#define LAPACK_EXPORT_slaqz0 F77_FUNC( slaqz0 , SLAQZ0 )
#define LAPACK_EXPORT_slaqz1 F77_FUNC( slaqz1 , SLAQZ1 )
#define LAPACK_EXPORT_slaqz2 F77_FUNC( slaqz2 , SLAQZ2 )
#define LAPACK_EXPORT_slaqz3 F77_FUNC( slaqz3 , SLAQZ3 )
#define LAPACK_EXPORT_slaqz4 F77_FUNC( slaqz4 , SLAQZ4 )
#define LAPACK_EXPORT_slarfb_gett F77_FUNC( slarfb_gett , SLARFB_GETT )
#define LAPACK_EXPORT_sorgtsqr_row F77_FUNC( sorgtsqr_row , SORGTSQR_ROW )
#define LAPACK_EXPORT_zgetsqrhrt F77_FUNC( zgetsqrhrt , ZGETSQRHRT )
#define LAPACK_EXPORT_zlaqz0 F77_FUNC( zlaqz0 , ZLAQZ0 )
#define LAPACK_EXPORT_zlaqz1 F77_FUNC( zlaqz1 , ZLAQZ1 )
#define LAPACK_EXPORT_zlaqz2 F77_FUNC( zlaqz2 , ZLAQZ2 )
#define LAPACK_EXPORT_zlaqz3 F77_FUNC( zlaqz3 , ZLAQZ3 )
#define LAPACK_EXPORT_zlarfb_gett F77_FUNC( zlarfb_gett , ZLARFB_GETT )
#define LAPACK_EXPORT_zungtsqr_row F77_FUNC( zungtsqr_row , ZUNGTSQR_ROW )
#define LAPACK_EXPORT_cbbcsd F77_FUNC( cbbcsd , CBBCSD )
#define LAPACK_EXPORT_cbdsqr F77_FUNC( cbdsqr , CBDSQR )
#define LAPACK_EXPORT_cgbbrd F77_FUNC( cgbbrd , CGBBRD )
#define LAPACK_EXPORT_cgbcon F77_FUNC( cgbcon , CGBCON )
#define LAPACK_EXPORT_cgbequb F77_FUNC( cgbequb , CGBEQUB )
#define LAPACK_EXPORT_cgbequ F77_FUNC( cgbequ , CGBEQU )
#define LAPACK_EXPORT_cgbrfs F77_FUNC( cgbrfs , CGBRFS )
#define LAPACK_EXPORT_cgbrfsx F77_FUNC( cgbrfsx , CGBRFSX )
#define LAPACK_EXPORT_cgbsv F77_FUNC( cgbsv , CGBSV )
#define LAPACK_EXPORT_cgbsvx F77_FUNC( cgbsvx , CGBSVX )
#define LAPACK_EXPORT_cgbsvxx F77_FUNC( cgbsvxx , CGBSVXX )
#define LAPACK_EXPORT_cgbtf2 F77_FUNC( cgbtf2 , CGBTF2 )
#define LAPACK_EXPORT_cgbtrf F77_FUNC( cgbtrf , CGBTRF )
#define LAPACK_EXPORT_cgbtrs F77_FUNC( cgbtrs , CGBTRS )
#define LAPACK_EXPORT_cgebak F77_FUNC( cgebak , CGEBAK )
#define LAPACK_EXPORT_cgebal F77_FUNC( cgebal , CGEBAL )
#define LAPACK_EXPORT_cgebd2 F77_FUNC( cgebd2 , CGEBD2 )
#define LAPACK_EXPORT_cgebrd F77_FUNC( cgebrd , CGEBRD )
#define LAPACK_EXPORT_cgecon F77_FUNC( cgecon , CGECON )
#define LAPACK_EXPORT_cgeequb F77_FUNC( cgeequb , CGEEQUB )
#define LAPACK_EXPORT_cgeequ F77_FUNC( cgeequ , CGEEQU )
#define LAPACK_EXPORT_cgees F77_FUNC( cgees , CGEES )
#define LAPACK_EXPORT_cgeesx F77_FUNC( cgeesx , CGEESX )
#define LAPACK_EXPORT_cgeev F77_FUNC( cgeev , CGEEV )
#define LAPACK_EXPORT_cgeevx F77_FUNC( cgeevx , CGEEVX )
#define LAPACK_EXPORT_cgegs F77_FUNC( cgegs , CGEGS )
#define LAPACK_EXPORT_cgegv F77_FUNC( cgegv , CGEGV )
#define LAPACK_EXPORT_cgehd2 F77_FUNC( cgehd2 , CGEHD2 )
#define LAPACK_EXPORT_cgehrd F77_FUNC( cgehrd , CGEHRD )
#define LAPACK_EXPORT_cgejsv F77_FUNC( cgejsv , CGEJSV )
#define LAPACK_EXPORT_cgelq2 F77_FUNC( cgelq2 , CGELQ2 )
#define LAPACK_EXPORT_cgelq F77_FUNC( cgelq , CGELQ )
#define LAPACK_EXPORT_cgelqf F77_FUNC( cgelqf , CGELQF )
#define LAPACK_EXPORT_cgelqt3 F77_FUNC( cgelqt3 , CGELQT3 )
#define LAPACK_EXPORT_cgelqt F77_FUNC( cgelqt , CGELQT )
#define LAPACK_EXPORT_cgelsd F77_FUNC( cgelsd , CGELSD )
#define LAPACK_EXPORT_cgels F77_FUNC( cgels , CGELS )
#define LAPACK_EXPORT_cgelss F77_FUNC( cgelss , CGELSS )
#define LAPACK_EXPORT_cgelsx F77_FUNC( cgelsx , CGELSX )
#define LAPACK_EXPORT_cgelsy F77_FUNC( cgelsy , CGELSY )
#define LAPACK_EXPORT_cgemlq F77_FUNC( cgemlq , CGEMLQ )
#define LAPACK_EXPORT_cgemlqt F77_FUNC( cgemlqt , CGEMLQT )
#define LAPACK_EXPORT_cgemqr F77_FUNC( cgemqr , CGEMQR )
#define LAPACK_EXPORT_cgemqrt F77_FUNC( cgemqrt , CGEMQRT )
#define LAPACK_EXPORT_cgeql2 F77_FUNC( cgeql2 , CGEQL2 )
#define LAPACK_EXPORT_cgeqlf F77_FUNC( cgeqlf , CGEQLF )
#define LAPACK_EXPORT_cgeqp3 F77_FUNC( cgeqp3 , CGEQP3 )
#define LAPACK_EXPORT_cgeqpf F77_FUNC( cgeqpf , CGEQPF )
#define LAPACK_EXPORT_cgeqr2 F77_FUNC( cgeqr2 , CGEQR2 )
#define LAPACK_EXPORT_cgeqr2p F77_FUNC( cgeqr2p , CGEQR2P )
#define LAPACK_EXPORT_cgeqr F77_FUNC( cgeqr , CGEQR )
#define LAPACK_EXPORT_cgeqrf F77_FUNC( cgeqrf , CGEQRF )
#define LAPACK_EXPORT_cgeqrfp F77_FUNC( cgeqrfp , CGEQRFP )
#define LAPACK_EXPORT_cgeqrt2 F77_FUNC( cgeqrt2 , CGEQRT2 )
#define LAPACK_EXPORT_cgeqrt3 F77_FUNC( cgeqrt3 , CGEQRT3 )
#define LAPACK_EXPORT_cgeqrt F77_FUNC( cgeqrt , CGEQRT )
#define LAPACK_EXPORT_cgerfs F77_FUNC( cgerfs , CGERFS )
#define LAPACK_EXPORT_cgerfsx F77_FUNC( cgerfsx , CGERFSX )
#define LAPACK_EXPORT_cgerq2 F77_FUNC( cgerq2 , CGERQ2 )
#define LAPACK_EXPORT_cgerqf F77_FUNC( cgerqf , CGERQF )
#define LAPACK_EXPORT_cgesc2 F77_FUNC( cgesc2 , CGESC2 )
#define LAPACK_EXPORT_cgesdd F77_FUNC( cgesdd , CGESDD )
#define LAPACK_EXPORT_cgesvd F77_FUNC( cgesvd , CGESVD )
#define LAPACK_EXPORT_cgesvdq F77_FUNC( cgesvdq , CGESVDQ )
#define LAPACK_EXPORT_cgesvdx F77_FUNC( cgesvdx , CGESVDX )
#define LAPACK_EXPORT_cgesv F77_FUNC( cgesv , CGESV )
#define LAPACK_EXPORT_cgesvj F77_FUNC( cgesvj , CGESVJ )
#define LAPACK_EXPORT_cgesvx F77_FUNC( cgesvx , CGESVX )
#define LAPACK_EXPORT_cgesvxx F77_FUNC( cgesvxx , CGESVXX )
#define LAPACK_EXPORT_cgetc2 F77_FUNC( cgetc2 , CGETC2 )
#define LAPACK_EXPORT_cgetf2 F77_FUNC( cgetf2 , CGETF2 )
#define LAPACK_EXPORT_cgetrf2 F77_FUNC( cgetrf2 , CGETRF2 )
#define LAPACK_EXPORT_cgetrf F77_FUNC( cgetrf , CGETRF )
#define LAPACK_EXPORT_cgetri F77_FUNC( cgetri , CGETRI )
#define LAPACK_EXPORT_cgetrs F77_FUNC( cgetrs , CGETRS )
#define LAPACK_EXPORT_cgetsls F77_FUNC( cgetsls , CGETSLS )
#define LAPACK_EXPORT_cggbak F77_FUNC( cggbak , CGGBAK )
#define LAPACK_EXPORT_cggbal F77_FUNC( cggbal , CGGBAL )
#define LAPACK_EXPORT_cgges3 F77_FUNC( cgges3 , CGGES3 )
#define LAPACK_EXPORT_cgges F77_FUNC( cgges , CGGES )
#define LAPACK_EXPORT_cggesx F77_FUNC( cggesx , CGGESX )
#define LAPACK_EXPORT_cggev3 F77_FUNC( cggev3 , CGGEV3 )
#define LAPACK_EXPORT_cggev F77_FUNC( cggev , CGGEV )
#define LAPACK_EXPORT_cggevx F77_FUNC( cggevx , CGGEVX )
#define LAPACK_EXPORT_cggglm F77_FUNC( cggglm , CGGGLM )
#define LAPACK_EXPORT_cgghd3 F77_FUNC( cgghd3 , CGGHD3 )
#define LAPACK_EXPORT_cgghrd F77_FUNC( cgghrd , CGGHRD )
#define LAPACK_EXPORT_cgglse F77_FUNC( cgglse , CGGLSE )
#define LAPACK_EXPORT_cggqrf F77_FUNC( cggqrf , CGGQRF )
#define LAPACK_EXPORT_cggrqf F77_FUNC( cggrqf , CGGRQF )
#define LAPACK_EXPORT_cggsvd3 F77_FUNC( cggsvd3 , CGGSVD3 )
#define LAPACK_EXPORT_cggsvd F77_FUNC( cggsvd , CGGSVD )
#define LAPACK_EXPORT_cggsvp3 F77_FUNC( cggsvp3 , CGGSVP3 )
#define LAPACK_EXPORT_cggsvp F77_FUNC( cggsvp , CGGSVP )
#define LAPACK_EXPORT_cgsvj0 F77_FUNC( cgsvj0 , CGSVJ0 )
#define LAPACK_EXPORT_cgsvj1 F77_FUNC( cgsvj1 , CGSVJ1 )
#define LAPACK_EXPORT_cgtcon F77_FUNC( cgtcon , CGTCON )
#define LAPACK_EXPORT_cgtrfs F77_FUNC( cgtrfs , CGTRFS )
#define LAPACK_EXPORT_cgtsv F77_FUNC( cgtsv , CGTSV )
#define LAPACK_EXPORT_cgtsvx F77_FUNC( cgtsvx , CGTSVX )
#define LAPACK_EXPORT_cgttrf F77_FUNC( cgttrf , CGTTRF )
#define LAPACK_EXPORT_cgttrs F77_FUNC( cgttrs , CGTTRS )
#define LAPACK_EXPORT_cgtts2 F77_FUNC( cgtts2 , CGTTS2 )
#define LAPACK_EXPORT_chb2st_kernels F77_FUNC( chb2st_kernels , CHB2ST_KERNELS )
#define LAPACK_EXPORT_chbev_2stage F77_FUNC( chbev_2stage , CHBEV_2STAGE )
#define LAPACK_EXPORT_chbevd_2stage F77_FUNC( chbevd_2stage , CHBEVD_2STAGE )
#define LAPACK_EXPORT_chbevd F77_FUNC( chbevd , CHBEVD )
#define LAPACK_EXPORT_chbev F77_FUNC( chbev , CHBEV )
#define LAPACK_EXPORT_chbevx_2stage F77_FUNC( chbevx_2stage , CHBEVX_2STAGE )
#define LAPACK_EXPORT_chbevx F77_FUNC( chbevx , CHBEVX )
#define LAPACK_EXPORT_chbgst F77_FUNC( chbgst , CHBGST )
#define LAPACK_EXPORT_chbgvd F77_FUNC( chbgvd , CHBGVD )
#define LAPACK_EXPORT_chbgv F77_FUNC( chbgv , CHBGV )
#define LAPACK_EXPORT_chbgvx F77_FUNC( chbgvx , CHBGVX )
#define LAPACK_EXPORT_chbtrd F77_FUNC( chbtrd , CHBTRD )
#define LAPACK_EXPORT_checon_3 F77_FUNC( checon_3 , CHECON_3 )
#define LAPACK_EXPORT_checon F77_FUNC( checon , CHECON )
#define LAPACK_EXPORT_checon_rook F77_FUNC( checon_rook , CHECON_ROOK )
#define LAPACK_EXPORT_cheequb F77_FUNC( cheequb , CHEEQUB )
#define LAPACK_EXPORT_cheev_2stage F77_FUNC( cheev_2stage , CHEEV_2STAGE )
#define LAPACK_EXPORT_cheevd_2stage F77_FUNC( cheevd_2stage , CHEEVD_2STAGE )
#define LAPACK_EXPORT_cheevd F77_FUNC( cheevd , CHEEVD )
#define LAPACK_EXPORT_cheev  F77_FUNC( cheev  , CHEEV  )
#define LAPACK_EXPORT_cheev F77_FUNC( cheev , CHEEV )
#define LAPACK_EXPORT_cheevr_2stage F77_FUNC( cheevr_2stage , CHEEVR_2STAGE )
#define LAPACK_EXPORT_cheevr F77_FUNC( cheevr , CHEEVR )
#define LAPACK_EXPORT_cheevx_2stage F77_FUNC( cheevx_2stage , CHEEVX_2STAGE )
#define LAPACK_EXPORT_cheevx F77_FUNC( cheevx , CHEEVX )
#define LAPACK_EXPORT_chegs2 F77_FUNC( chegs2 , CHEGS2 )
#define LAPACK_EXPORT_chegst F77_FUNC( chegst , CHEGST )
#define LAPACK_EXPORT_chegv_2stage F77_FUNC( chegv_2stage , CHEGV_2STAGE )
#define LAPACK_EXPORT_chegvd F77_FUNC( chegvd , CHEGVD )
#define LAPACK_EXPORT_chegv F77_FUNC( chegv , CHEGV )
#define LAPACK_EXPORT_chegvx F77_FUNC( chegvx , CHEGVX )
#define LAPACK_EXPORT_cherfs F77_FUNC( cherfs , CHERFS )
#define LAPACK_EXPORT_cherfsx F77_FUNC( cherfsx , CHERFSX )
#define LAPACK_EXPORT_chesv_aa_2stage F77_FUNC( chesv_aa_2stage , CHESV_AA_2STAGE )
#define LAPACK_EXPORT_chesv_aa F77_FUNC( chesv_aa , CHESV_AA )
#define LAPACK_EXPORT_chesv F77_FUNC( chesv , CHESV )
#define LAPACK_EXPORT_chesv_rk F77_FUNC( chesv_rk , CHESV_RK )
#define LAPACK_EXPORT_chesv_rook F77_FUNC( chesv_rook , CHESV_ROOK )
#define LAPACK_EXPORT_chesvx F77_FUNC( chesvx , CHESVX )
#define LAPACK_EXPORT_chesvxx F77_FUNC( chesvxx , CHESVXX )
#define LAPACK_EXPORT_cheswapr F77_FUNC( cheswapr , CHESWAPR )
#define LAPACK_EXPORT_chetd2 F77_FUNC( chetd2 , CHETD2 )
#define LAPACK_EXPORT_chetf2 F77_FUNC( chetf2 , CHETF2 )
#define LAPACK_EXPORT_chetf2_rk F77_FUNC( chetf2_rk , CHETF2_RK )
#define LAPACK_EXPORT_chetf2_rook F77_FUNC( chetf2_rook , CHETF2_ROOK )
#define LAPACK_EXPORT_chetrd_2stage F77_FUNC( chetrd_2stage , CHETRD_2STAGE )
#define LAPACK_EXPORT_chetrd F77_FUNC( chetrd , CHETRD )
#define LAPACK_EXPORT_chetrd_hb2st F77_FUNC( chetrd_hb2st , CHETRD_HB2ST )
#define LAPACK_EXPORT_chetrd_he2hb F77_FUNC( chetrd_he2hb , CHETRD_HE2HB )
#define LAPACK_EXPORT_chetrf_aa_2stage F77_FUNC( chetrf_aa_2stage , CHETRF_AA_2STAGE )
#define LAPACK_EXPORT_chetrf_aa F77_FUNC( chetrf_aa , CHETRF_AA )
#define LAPACK_EXPORT_chetrf F77_FUNC( chetrf , CHETRF )
#define LAPACK_EXPORT_chetrf_rk F77_FUNC( chetrf_rk , CHETRF_RK )
#define LAPACK_EXPORT_chetrf_rook F77_FUNC( chetrf_rook , CHETRF_ROOK )
#define LAPACK_EXPORT_chetri2 F77_FUNC( chetri2 , CHETRI2 )
#define LAPACK_EXPORT_chetri2x F77_FUNC( chetri2x , CHETRI2X )
#define LAPACK_EXPORT_chetri_3 F77_FUNC( chetri_3 , CHETRI_3 )
#define LAPACK_EXPORT_chetri_3x F77_FUNC( chetri_3x , CHETRI_3X )
#define LAPACK_EXPORT_chetri F77_FUNC( chetri , CHETRI )
#define LAPACK_EXPORT_chetri_rook F77_FUNC( chetri_rook , CHETRI_ROOK )
#define LAPACK_EXPORT_chetrs2 F77_FUNC( chetrs2 , CHETRS2 )
#define LAPACK_EXPORT_chetrs_3 F77_FUNC( chetrs_3 , CHETRS_3 )
#define LAPACK_EXPORT_chetrs_aa_2stage F77_FUNC( chetrs_aa_2stage , CHETRS_AA_2STAGE )
#define LAPACK_EXPORT_chetrs_aa F77_FUNC( chetrs_aa , CHETRS_AA )
#define LAPACK_EXPORT_chetrs F77_FUNC( chetrs , CHETRS )
#define LAPACK_EXPORT_chetrs_rook F77_FUNC( chetrs_rook , CHETRS_ROOK )
#define LAPACK_EXPORT_chfrk F77_FUNC( chfrk , CHFRK )
#define LAPACK_EXPORT_chgeqz F77_FUNC( chgeqz , CHGEQZ )
#define LAPACK_EXPORT_chla_transtype F77_FUNC( chla_transtype , CHLA_TRANSTYPE )
#define LAPACK_EXPORT_chpcon F77_FUNC( chpcon , CHPCON )
#define LAPACK_EXPORT_chpevd F77_FUNC( chpevd , CHPEVD )
#define LAPACK_EXPORT_chpev F77_FUNC( chpev , CHPEV )
#define LAPACK_EXPORT_chpevx F77_FUNC( chpevx , CHPEVX )
#define LAPACK_EXPORT_chpgst F77_FUNC( chpgst , CHPGST )
#define LAPACK_EXPORT_chpgvd F77_FUNC( chpgvd , CHPGVD )
#define LAPACK_EXPORT_chpgv F77_FUNC( chpgv , CHPGV )
#define LAPACK_EXPORT_chpgvx F77_FUNC( chpgvx , CHPGVX )
#define LAPACK_EXPORT_chprfs F77_FUNC( chprfs , CHPRFS )
#define LAPACK_EXPORT_chpsv F77_FUNC( chpsv , CHPSV )
#define LAPACK_EXPORT_chpsvx F77_FUNC( chpsvx , CHPSVX )
#define LAPACK_EXPORT_chptrd F77_FUNC( chptrd , CHPTRD )
#define LAPACK_EXPORT_chptrf F77_FUNC( chptrf , CHPTRF )
#define LAPACK_EXPORT_chptri F77_FUNC( chptri , CHPTRI )
#define LAPACK_EXPORT_chptrs F77_FUNC( chptrs , CHPTRS )
#define LAPACK_EXPORT_chsein F77_FUNC( chsein , CHSEIN )
#define LAPACK_EXPORT_chseqr F77_FUNC( chseqr , CHSEQR )
#define LAPACK_EXPORT_clabrd F77_FUNC( clabrd , CLABRD )
#define LAPACK_EXPORT_clacgv F77_FUNC( clacgv , CLACGV )
#define LAPACK_EXPORT_clacn2 F77_FUNC( clacn2 , CLACN2 )
#define LAPACK_EXPORT_clacon F77_FUNC( clacon , CLACON )
#define LAPACK_EXPORT_clacp2 F77_FUNC( clacp2 , CLACP2 )
#define LAPACK_EXPORT_clacpy F77_FUNC( clacpy , CLACPY )
#define LAPACK_EXPORT_clacrm F77_FUNC( clacrm , CLACRM )
#define LAPACK_EXPORT_clacrt F77_FUNC( clacrt , CLACRT )
#define LAPACK_EXPORT_cladiv F77_FUNC( cladiv , CLADIV )
#define LAPACK_EXPORT_claed0 F77_FUNC( claed0 , CLAED0 )
#define LAPACK_EXPORT_claed7 F77_FUNC( claed7 , CLAED7 )
#define LAPACK_EXPORT_claed8 F77_FUNC( claed8 , CLAED8 )
#define LAPACK_EXPORT_claein F77_FUNC( claein , CLAEIN )
#define LAPACK_EXPORT_claesy F77_FUNC( claesy , CLAESY )
#define LAPACK_EXPORT_claev2 F77_FUNC( claev2 , CLAEV2 )
#define LAPACK_EXPORT_clag2z F77_FUNC( clag2z , CLAG2Z )
#define LAPACK_EXPORT_cla_gbamv F77_FUNC( cla_gbamv , CLA_GBAMV )
#define LAPACK_EXPORT_cla_gbrcond_c F77_FUNC( cla_gbrcond_c , CLA_GBRCOND_C )
#define LAPACK_EXPORT_cla_gbrcond_x F77_FUNC( cla_gbrcond_x , CLA_GBRCOND_X )
#define LAPACK_EXPORT_cla_gbrfsx_extended F77_FUNC( cla_gbrfsx_extended , CLA_GBRFSX_EXTENDED )
#define LAPACK_EXPORT_cla_gbrpvgrw F77_FUNC( cla_gbrpvgrw , CLA_GBRPVGRW )
#define LAPACK_EXPORT_cla_geamv F77_FUNC( cla_geamv , CLA_GEAMV )
#define LAPACK_EXPORT_cla_gercond_c F77_FUNC( cla_gercond_c , CLA_GERCOND_C )
#define LAPACK_EXPORT_cla_gercond_x F77_FUNC( cla_gercond_x , CLA_GERCOND_X )
#define LAPACK_EXPORT_cla_gerfsx_extended F77_FUNC( cla_gerfsx_extended , CLA_GERFSX_EXTENDED )
#define LAPACK_EXPORT_cla_gerpvgrw F77_FUNC( cla_gerpvgrw , CLA_GERPVGRW )
#define LAPACK_EXPORT_clags2 F77_FUNC( clags2 , CLAGS2 )
#define LAPACK_EXPORT_clagtm F77_FUNC( clagtm , CLAGTM )
#define LAPACK_EXPORT_cla_heamv F77_FUNC( cla_heamv , CLA_HEAMV )
#define LAPACK_EXPORT_clahef_aa F77_FUNC( clahef_aa , CLAHEF_AA )
#define LAPACK_EXPORT_clahef F77_FUNC( clahef , CLAHEF )
#define LAPACK_EXPORT_clahef_rk F77_FUNC( clahef_rk , CLAHEF_RK )
#define LAPACK_EXPORT_clahef_rook F77_FUNC( clahef_rook , CLAHEF_ROOK )
#define LAPACK_EXPORT_cla_hercond_c F77_FUNC( cla_hercond_c , CLA_HERCOND_C )
#define LAPACK_EXPORT_cla_hercond_x F77_FUNC( cla_hercond_x , CLA_HERCOND_X )
#define LAPACK_EXPORT_cla_herfsx_extended F77_FUNC( cla_herfsx_extended , CLA_HERFSX_EXTENDED )
#define LAPACK_EXPORT_cla_herpvgrw F77_FUNC( cla_herpvgrw , CLA_HERPVGRW )
#define LAPACK_EXPORT_clahqr F77_FUNC( clahqr , CLAHQR )
#define LAPACK_EXPORT_clahr2 F77_FUNC( clahr2 , CLAHR2 )
#define LAPACK_EXPORT_clahrd F77_FUNC( clahrd , CLAHRD )
#define LAPACK_EXPORT_claic1 F77_FUNC( claic1 , CLAIC1 )
#define LAPACK_EXPORT_cla_lin_berr F77_FUNC( cla_lin_berr , CLA_LIN_BERR )
#define LAPACK_EXPORT_clals0 F77_FUNC( clals0 , CLALS0 )
#define LAPACK_EXPORT_clalsa F77_FUNC( clalsa , CLALSA )
#define LAPACK_EXPORT_clalsd F77_FUNC( clalsd , CLALSD )
#define LAPACK_EXPORT_clamswlq F77_FUNC( clamswlq , CLAMSWLQ )
#define LAPACK_EXPORT_clamtsqr F77_FUNC( clamtsqr , CLAMTSQR )
#define LAPACK_EXPORT_clangb F77_FUNC( clangb , CLANGB )
#define LAPACK_EXPORT_clange F77_FUNC( clange , CLANGE )
#define LAPACK_EXPORT_clangt F77_FUNC( clangt , CLANGT )
#define LAPACK_EXPORT_clanhb F77_FUNC( clanhb , CLANHB )
#define LAPACK_EXPORT_clanhe F77_FUNC( clanhe , CLANHE )
#define LAPACK_EXPORT_clanhf F77_FUNC( clanhf , CLANHF )
#define LAPACK_EXPORT_clanhp F77_FUNC( clanhp , CLANHP )
#define LAPACK_EXPORT_clanhs F77_FUNC( clanhs , CLANHS )
#define LAPACK_EXPORT_clanht F77_FUNC( clanht , CLANHT )
#define LAPACK_EXPORT_clansb F77_FUNC( clansb , CLANSB )
#define LAPACK_EXPORT_clansp F77_FUNC( clansp , CLANSP )
#define LAPACK_EXPORT_clansy F77_FUNC( clansy , CLANSY )
#define LAPACK_EXPORT_clantb F77_FUNC( clantb , CLANTB )
#define LAPACK_EXPORT_clantp F77_FUNC( clantp , CLANTP )
#define LAPACK_EXPORT_clantr F77_FUNC( clantr , CLANTR )
#define LAPACK_EXPORT_clapll F77_FUNC( clapll , CLAPLL )
#define LAPACK_EXPORT_clapmr F77_FUNC( clapmr , CLAPMR )
#define LAPACK_EXPORT_clapmt F77_FUNC( clapmt , CLAPMT )
#define LAPACK_EXPORT_cla_porcond_c F77_FUNC( cla_porcond_c , CLA_PORCOND_C )
#define LAPACK_EXPORT_cla_porcond_x F77_FUNC( cla_porcond_x , CLA_PORCOND_X )
#define LAPACK_EXPORT_cla_porfsx_extended F77_FUNC( cla_porfsx_extended , CLA_PORFSX_EXTENDED )
#define LAPACK_EXPORT_cla_porpvgrw F77_FUNC( cla_porpvgrw , CLA_PORPVGRW )
#define LAPACK_EXPORT_claqgb F77_FUNC( claqgb , CLAQGB )
#define LAPACK_EXPORT_claqge F77_FUNC( claqge , CLAQGE )
#define LAPACK_EXPORT_claqhb F77_FUNC( claqhb , CLAQHB )
#define LAPACK_EXPORT_claqhe F77_FUNC( claqhe , CLAQHE )
#define LAPACK_EXPORT_claqhp F77_FUNC( claqhp , CLAQHP )
#define LAPACK_EXPORT_claqp2 F77_FUNC( claqp2 , CLAQP2 )
#define LAPACK_EXPORT_claqps F77_FUNC( claqps , CLAQPS )
#define LAPACK_EXPORT_claqr0 F77_FUNC( claqr0 , CLAQR0 )
#define LAPACK_EXPORT_claqr1 F77_FUNC( claqr1 , CLAQR1 )
#define LAPACK_EXPORT_claqr2 F77_FUNC( claqr2 , CLAQR2 )
#define LAPACK_EXPORT_claqr3 F77_FUNC( claqr3 , CLAQR3 )
#define LAPACK_EXPORT_claqr4 F77_FUNC( claqr4 , CLAQR4 )
#define LAPACK_EXPORT_claqr5 F77_FUNC( claqr5 , CLAQR5 )
#define LAPACK_EXPORT_claqsb F77_FUNC( claqsb , CLAQSB )
#define LAPACK_EXPORT_claqsp F77_FUNC( claqsp , CLAQSP )
#define LAPACK_EXPORT_claqsy F77_FUNC( claqsy , CLAQSY )
#define LAPACK_EXPORT_clar1v F77_FUNC( clar1v , CLAR1V )
#define LAPACK_EXPORT_clar2v F77_FUNC( clar2v , CLAR2V )
#define LAPACK_EXPORT_clarcm F77_FUNC( clarcm , CLARCM )
#define LAPACK_EXPORT_clarfb F77_FUNC( clarfb , CLARFB )
#define LAPACK_EXPORT_clarf F77_FUNC( clarf , CLARF )
#define LAPACK_EXPORT_clarfg F77_FUNC( clarfg , CLARFG )
#define LAPACK_EXPORT_clarfgp F77_FUNC( clarfgp , CLARFGP )
#define LAPACK_EXPORT_clarft F77_FUNC( clarft , CLARFT )
#define LAPACK_EXPORT_clarfx F77_FUNC( clarfx , CLARFX )
#define LAPACK_EXPORT_clarfy F77_FUNC( clarfy , CLARFY )
#define LAPACK_EXPORT_clargv F77_FUNC( clargv , CLARGV )
#define LAPACK_EXPORT_clarnv F77_FUNC( clarnv , CLARNV )
#define LAPACK_EXPORT_clarrv F77_FUNC( clarrv , CLARRV )
#define LAPACK_EXPORT_clarscl2 F77_FUNC( clarscl2 , CLARSCL2 )
#define LAPACK_EXPORT_clartg F77_FUNC( clartg , CLARTG )
#define LAPACK_EXPORT_clartv F77_FUNC( clartv , CLARTV )
#define LAPACK_EXPORT_clarzb F77_FUNC( clarzb , CLARZB )
#define LAPACK_EXPORT_clarz F77_FUNC( clarz , CLARZ )
#define LAPACK_EXPORT_clarzt F77_FUNC( clarzt , CLARZT )
#define LAPACK_EXPORT_clascl2 F77_FUNC( clascl2 , CLASCL2 )
#define LAPACK_EXPORT_clascl F77_FUNC( clascl , CLASCL )
#define LAPACK_EXPORT_claset F77_FUNC( claset , CLASET )
#define LAPACK_EXPORT_clasr F77_FUNC( clasr , CLASR )
#define LAPACK_EXPORT_classq F77_FUNC( classq , CLASSQ )
#define LAPACK_EXPORT_claswlq F77_FUNC( claswlq , CLASWLQ )
#define LAPACK_EXPORT_claswp F77_FUNC( claswp , CLASWP )
#define LAPACK_EXPORT_cla_syamv F77_FUNC( cla_syamv , CLA_SYAMV )
#define LAPACK_EXPORT_clasyf_aa F77_FUNC( clasyf_aa , CLASYF_AA )
#define LAPACK_EXPORT_clasyf F77_FUNC( clasyf , CLASYF )
#define LAPACK_EXPORT_clasyf_rk F77_FUNC( clasyf_rk , CLASYF_RK )
#define LAPACK_EXPORT_clasyf_rook F77_FUNC( clasyf_rook , CLASYF_ROOK )
#define LAPACK_EXPORT_cla_syrcond_c F77_FUNC( cla_syrcond_c , CLA_SYRCOND_C )
#define LAPACK_EXPORT_cla_syrcond_x F77_FUNC( cla_syrcond_x , CLA_SYRCOND_X )
#define LAPACK_EXPORT_cla_syrfsx_extended F77_FUNC( cla_syrfsx_extended , CLA_SYRFSX_EXTENDED )
#define LAPACK_EXPORT_cla_syrpvgrw F77_FUNC( cla_syrpvgrw , CLA_SYRPVGRW )
#define LAPACK_EXPORT_clatbs F77_FUNC( clatbs , CLATBS )
#define LAPACK_EXPORT_clatdf F77_FUNC( clatdf , CLATDF )
#define LAPACK_EXPORT_clatps F77_FUNC( clatps , CLATPS )
#define LAPACK_EXPORT_clatrd F77_FUNC( clatrd , CLATRD )
#define LAPACK_EXPORT_clatrs F77_FUNC( clatrs , CLATRS )
#define LAPACK_EXPORT_clatrz F77_FUNC( clatrz , CLATRZ )
#define LAPACK_EXPORT_clatsqr F77_FUNC( clatsqr , CLATSQR )
#define LAPACK_EXPORT_clatzm F77_FUNC( clatzm , CLATZM )
#define LAPACK_EXPORT_claunhr_col_getrfnp2 F77_FUNC( claunhr_col_getrfnp2 , CLAUNHR_COL_GETRFNP2 )
#define LAPACK_EXPORT_claunhr_col_getrfnp F77_FUNC( claunhr_col_getrfnp , CLAUNHR_COL_GETRFNP )
#define LAPACK_EXPORT_clauu2 F77_FUNC( clauu2 , CLAUU2 )
#define LAPACK_EXPORT_clauum F77_FUNC( clauum , CLAUUM )
#define LAPACK_EXPORT_cla_wwaddw F77_FUNC( cla_wwaddw , CLA_WWADDW )
#define LAPACK_EXPORT_cpbcon F77_FUNC( cpbcon , CPBCON )
#define LAPACK_EXPORT_cpbequ F77_FUNC( cpbequ , CPBEQU )
#define LAPACK_EXPORT_cpbrfs F77_FUNC( cpbrfs , CPBRFS )
#define LAPACK_EXPORT_cpbstf F77_FUNC( cpbstf , CPBSTF )
#define LAPACK_EXPORT_cpbsv F77_FUNC( cpbsv , CPBSV )
#define LAPACK_EXPORT_cpbsvx F77_FUNC( cpbsvx , CPBSVX )
#define LAPACK_EXPORT_cpbtf2 F77_FUNC( cpbtf2 , CPBTF2 )
#define LAPACK_EXPORT_cpbtrf F77_FUNC( cpbtrf , CPBTRF )
#define LAPACK_EXPORT_cpbtrs F77_FUNC( cpbtrs , CPBTRS )
#define LAPACK_EXPORT_cpftrf F77_FUNC( cpftrf , CPFTRF )
#define LAPACK_EXPORT_cpftri F77_FUNC( cpftri , CPFTRI )
#define LAPACK_EXPORT_cpftrs F77_FUNC( cpftrs , CPFTRS )
#define LAPACK_EXPORT_cpocon F77_FUNC( cpocon , CPOCON )
#define LAPACK_EXPORT_cpoequb F77_FUNC( cpoequb , CPOEQUB )
#define LAPACK_EXPORT_cpoequ F77_FUNC( cpoequ , CPOEQU )
#define LAPACK_EXPORT_cporfs F77_FUNC( cporfs , CPORFS )
#define LAPACK_EXPORT_cporfsx F77_FUNC( cporfsx , CPORFSX )
#define LAPACK_EXPORT_cposv F77_FUNC( cposv , CPOSV )
#define LAPACK_EXPORT_cposvx F77_FUNC( cposvx , CPOSVX )
#define LAPACK_EXPORT_cposvxx F77_FUNC( cposvxx , CPOSVXX )
#define LAPACK_EXPORT_cpotf2 F77_FUNC( cpotf2 , CPOTF2 )
#define LAPACK_EXPORT_cpotrf2 F77_FUNC( cpotrf2 , CPOTRF2 )
#define LAPACK_EXPORT_cpotrf F77_FUNC( cpotrf , CPOTRF )
#define LAPACK_EXPORT_cpotri F77_FUNC( cpotri , CPOTRI )
#define LAPACK_EXPORT_cpotrs F77_FUNC( cpotrs , CPOTRS )
#define LAPACK_EXPORT_cppcon F77_FUNC( cppcon , CPPCON )
#define LAPACK_EXPORT_cppequ F77_FUNC( cppequ , CPPEQU )
#define LAPACK_EXPORT_cpprfs F77_FUNC( cpprfs , CPPRFS )
#define LAPACK_EXPORT_cppsv F77_FUNC( cppsv , CPPSV )
#define LAPACK_EXPORT_cppsvx F77_FUNC( cppsvx , CPPSVX )
#define LAPACK_EXPORT_cpptrf F77_FUNC( cpptrf , CPPTRF )
#define LAPACK_EXPORT_cpptri F77_FUNC( cpptri , CPPTRI )
#define LAPACK_EXPORT_cpptrs F77_FUNC( cpptrs , CPPTRS )
#define LAPACK_EXPORT_cpstf2 F77_FUNC( cpstf2 , CPSTF2 )
#define LAPACK_EXPORT_cpstrf F77_FUNC( cpstrf , CPSTRF )
#define LAPACK_EXPORT_cptcon F77_FUNC( cptcon , CPTCON )
#define LAPACK_EXPORT_cpteqr F77_FUNC( cpteqr , CPTEQR )
#define LAPACK_EXPORT_cptrfs F77_FUNC( cptrfs , CPTRFS )
#define LAPACK_EXPORT_cptsv F77_FUNC( cptsv , CPTSV )
#define LAPACK_EXPORT_cptsvx F77_FUNC( cptsvx , CPTSVX )
#define LAPACK_EXPORT_cpttrf F77_FUNC( cpttrf , CPTTRF )
#define LAPACK_EXPORT_cpttrs F77_FUNC( cpttrs , CPTTRS )
#define LAPACK_EXPORT_cptts2 F77_FUNC( cptts2 , CPTTS2 )
#define LAPACK_EXPORT_crot F77_FUNC( crot , CROT )
#define LAPACK_EXPORT_cspcon F77_FUNC( cspcon , CSPCON )
#define LAPACK_EXPORT_cspmv F77_FUNC( cspmv , CSPMV )
#define LAPACK_EXPORT_cspr F77_FUNC( cspr , CSPR )
#define LAPACK_EXPORT_csprfs F77_FUNC( csprfs , CSPRFS )
#define LAPACK_EXPORT_cspsv F77_FUNC( cspsv , CSPSV )
#define LAPACK_EXPORT_cspsvx F77_FUNC( cspsvx , CSPSVX )
#define LAPACK_EXPORT_csptrf F77_FUNC( csptrf , CSPTRF )
#define LAPACK_EXPORT_csptri F77_FUNC( csptri , CSPTRI )
#define LAPACK_EXPORT_csptrs F77_FUNC( csptrs , CSPTRS )
#define LAPACK_EXPORT_csrscl F77_FUNC( csrscl , CSRSCL )
#define LAPACK_EXPORT_cstedc F77_FUNC( cstedc , CSTEDC )
#define LAPACK_EXPORT_cstegr F77_FUNC( cstegr , CSTEGR )
#define LAPACK_EXPORT_cstein F77_FUNC( cstein , CSTEIN )
#define LAPACK_EXPORT_cstemr F77_FUNC( cstemr , CSTEMR )
#define LAPACK_EXPORT_csteqr F77_FUNC( csteqr , CSTEQR )
#define LAPACK_EXPORT_csycon_3 F77_FUNC( csycon_3 , CSYCON_3 )
#define LAPACK_EXPORT_csycon F77_FUNC( csycon , CSYCON )
#define LAPACK_EXPORT_csycon_rook F77_FUNC( csycon_rook , CSYCON_ROOK )
#define LAPACK_EXPORT_csyconv F77_FUNC( csyconv , CSYCONV )
#define LAPACK_EXPORT_csyconvf F77_FUNC( csyconvf , CSYCONVF )
#define LAPACK_EXPORT_csyconvf_rook F77_FUNC( csyconvf_rook , CSYCONVF_ROOK )
#define LAPACK_EXPORT_csyequb F77_FUNC( csyequb , CSYEQUB )
#define LAPACK_EXPORT_csymv F77_FUNC( csymv , CSYMV )
#define LAPACK_EXPORT_csyr F77_FUNC( csyr , CSYR )
#define LAPACK_EXPORT_csyrfs F77_FUNC( csyrfs , CSYRFS )
#define LAPACK_EXPORT_csyrfsx F77_FUNC( csyrfsx , CSYRFSX )
#define LAPACK_EXPORT_csysv_aa_2stage F77_FUNC( csysv_aa_2stage , CSYSV_AA_2STAGE )
#define LAPACK_EXPORT_csysv_aa F77_FUNC( csysv_aa , CSYSV_AA )
#define LAPACK_EXPORT_csysv F77_FUNC( csysv , CSYSV )
#define LAPACK_EXPORT_csysv_rk F77_FUNC( csysv_rk , CSYSV_RK )
#define LAPACK_EXPORT_csysv_rook F77_FUNC( csysv_rook , CSYSV_ROOK )
#define LAPACK_EXPORT_csysvx F77_FUNC( csysvx , CSYSVX )
#define LAPACK_EXPORT_csysvxx F77_FUNC( csysvxx , CSYSVXX )
#define LAPACK_EXPORT_csyswapr F77_FUNC( csyswapr , CSYSWAPR )
#define LAPACK_EXPORT_csytf2 F77_FUNC( csytf2 , CSYTF2 )
#define LAPACK_EXPORT_csytf2_rk F77_FUNC( csytf2_rk , CSYTF2_RK )
#define LAPACK_EXPORT_csytf2_rook F77_FUNC( csytf2_rook , CSYTF2_ROOK )
#define LAPACK_EXPORT_csytrf_aa_2stage F77_FUNC( csytrf_aa_2stage , CSYTRF_AA_2STAGE )
#define LAPACK_EXPORT_csytrf_aa F77_FUNC( csytrf_aa , CSYTRF_AA )
#define LAPACK_EXPORT_csytrf F77_FUNC( csytrf , CSYTRF )
#define LAPACK_EXPORT_csytrf_rk F77_FUNC( csytrf_rk , CSYTRF_RK )
#define LAPACK_EXPORT_csytrf_rook F77_FUNC( csytrf_rook , CSYTRF_ROOK )
#define LAPACK_EXPORT_csytri2 F77_FUNC( csytri2 , CSYTRI2 )
#define LAPACK_EXPORT_csytri2x F77_FUNC( csytri2x , CSYTRI2X )
#define LAPACK_EXPORT_csytri_3 F77_FUNC( csytri_3 , CSYTRI_3 )
#define LAPACK_EXPORT_csytri_3x F77_FUNC( csytri_3x , CSYTRI_3X )
#define LAPACK_EXPORT_csytri F77_FUNC( csytri , CSYTRI )
#define LAPACK_EXPORT_csytri_rook F77_FUNC( csytri_rook , CSYTRI_ROOK )
#define LAPACK_EXPORT_csytrs2 F77_FUNC( csytrs2 , CSYTRS2 )
#define LAPACK_EXPORT_csytrs_3 F77_FUNC( csytrs_3 , CSYTRS_3 )
#define LAPACK_EXPORT_csytrs_aa_2stage F77_FUNC( csytrs_aa_2stage , CSYTRS_AA_2STAGE )
#define LAPACK_EXPORT_csytrs_aa F77_FUNC( csytrs_aa , CSYTRS_AA )
#define LAPACK_EXPORT_csytrs F77_FUNC( csytrs , CSYTRS )
#define LAPACK_EXPORT_csytrs_rook F77_FUNC( csytrs_rook , CSYTRS_ROOK )
#define LAPACK_EXPORT_ctbcon F77_FUNC( ctbcon , CTBCON )
#define LAPACK_EXPORT_ctbrfs F77_FUNC( ctbrfs , CTBRFS )
#define LAPACK_EXPORT_ctbtrs F77_FUNC( ctbtrs , CTBTRS )
#define LAPACK_EXPORT_ctfsm F77_FUNC( ctfsm , CTFSM )
#define LAPACK_EXPORT_ctftri F77_FUNC( ctftri , CTFTRI )
#define LAPACK_EXPORT_ctfttp F77_FUNC( ctfttp , CTFTTP )
#define LAPACK_EXPORT_ctfttr F77_FUNC( ctfttr , CTFTTR )
#define LAPACK_EXPORT_ctgevc F77_FUNC( ctgevc , CTGEVC )
#define LAPACK_EXPORT_ctgex2 F77_FUNC( ctgex2 , CTGEX2 )
#define LAPACK_EXPORT_ctgexc F77_FUNC( ctgexc , CTGEXC )
#define LAPACK_EXPORT_ctgsen F77_FUNC( ctgsen , CTGSEN )
#define LAPACK_EXPORT_ctgsja F77_FUNC( ctgsja , CTGSJA )
#define LAPACK_EXPORT_ctgsna F77_FUNC( ctgsna , CTGSNA )
#define LAPACK_EXPORT_ctgsy2 F77_FUNC( ctgsy2 , CTGSY2 )
#define LAPACK_EXPORT_ctgsyl F77_FUNC( ctgsyl , CTGSYL )
#define LAPACK_EXPORT_ctpcon F77_FUNC( ctpcon , CTPCON )
#define LAPACK_EXPORT_ctplqt2 F77_FUNC( ctplqt2 , CTPLQT2 )
#define LAPACK_EXPORT_ctplqt F77_FUNC( ctplqt , CTPLQT )
#define LAPACK_EXPORT_ctpmlqt F77_FUNC( ctpmlqt , CTPMLQT )
#define LAPACK_EXPORT_ctpmqrt F77_FUNC( ctpmqrt , CTPMQRT )
#define LAPACK_EXPORT_ctpqrt2 F77_FUNC( ctpqrt2 , CTPQRT2 )
#define LAPACK_EXPORT_ctpqrt F77_FUNC( ctpqrt , CTPQRT )
#define LAPACK_EXPORT_ctprfb F77_FUNC( ctprfb , CTPRFB )
#define LAPACK_EXPORT_ctprfs F77_FUNC( ctprfs , CTPRFS )
#define LAPACK_EXPORT_ctptri F77_FUNC( ctptri , CTPTRI )
#define LAPACK_EXPORT_ctptrs F77_FUNC( ctptrs , CTPTRS )
#define LAPACK_EXPORT_ctpttf F77_FUNC( ctpttf , CTPTTF )
#define LAPACK_EXPORT_ctpttr F77_FUNC( ctpttr , CTPTTR )
#define LAPACK_EXPORT_ctrcon F77_FUNC( ctrcon , CTRCON )
#define LAPACK_EXPORT_ctrevc3 F77_FUNC( ctrevc3 , CTREVC3 )
#define LAPACK_EXPORT_ctrevc F77_FUNC( ctrevc , CTREVC )
#define LAPACK_EXPORT_ctrexc F77_FUNC( ctrexc , CTREXC )
#define LAPACK_EXPORT_ctrrfs F77_FUNC( ctrrfs , CTRRFS )
#define LAPACK_EXPORT_ctrsen F77_FUNC( ctrsen , CTRSEN )
#define LAPACK_EXPORT_ctrsna F77_FUNC( ctrsna , CTRSNA )
#define LAPACK_EXPORT_ctrsyl F77_FUNC( ctrsyl , CTRSYL )
#define LAPACK_EXPORT_ctrti2 F77_FUNC( ctrti2 , CTRTI2 )
#define LAPACK_EXPORT_ctrtri F77_FUNC( ctrtri , CTRTRI )
#define LAPACK_EXPORT_ctrtrs F77_FUNC( ctrtrs , CTRTRS )
#define LAPACK_EXPORT_ctrttf F77_FUNC( ctrttf , CTRTTF )
#define LAPACK_EXPORT_ctrttp F77_FUNC( ctrttp , CTRTTP )
#define LAPACK_EXPORT_ctzrqf F77_FUNC( ctzrqf , CTZRQF )
#define LAPACK_EXPORT_ctzrzf F77_FUNC( ctzrzf , CTZRZF )
#define LAPACK_EXPORT_cunbdb1 F77_FUNC( cunbdb1 , CUNBDB1 )
#define LAPACK_EXPORT_cunbdb2 F77_FUNC( cunbdb2 , CUNBDB2 )
#define LAPACK_EXPORT_cunbdb3 F77_FUNC( cunbdb3 , CUNBDB3 )
#define LAPACK_EXPORT_cunbdb4 F77_FUNC( cunbdb4 , CUNBDB4 )
#define LAPACK_EXPORT_cunbdb5 F77_FUNC( cunbdb5 , CUNBDB5 )
#define LAPACK_EXPORT_cunbdb6 F77_FUNC( cunbdb6 , CUNBDB6 )
#define LAPACK_EXPORT_cunbdb F77_FUNC( cunbdb , CUNBDB )
#define LAPACK_EXPORT_cuncsd2by1 F77_FUNC( cuncsd2by1 , CUNCSD2BY1 )
#define LAPACK_EXPORT_cuncsd F77_FUNC( cuncsd , CUNCSD )
#define LAPACK_EXPORT_cung2l F77_FUNC( cung2l , CUNG2L )
#define LAPACK_EXPORT_cung2r F77_FUNC( cung2r , CUNG2R )
#define LAPACK_EXPORT_cungbr F77_FUNC( cungbr , CUNGBR )
#define LAPACK_EXPORT_cunghr F77_FUNC( cunghr , CUNGHR )
#define LAPACK_EXPORT_cungl2 F77_FUNC( cungl2 , CUNGL2 )
#define LAPACK_EXPORT_cunglq F77_FUNC( cunglq , CUNGLQ )
#define LAPACK_EXPORT_cungql F77_FUNC( cungql , CUNGQL )
#define LAPACK_EXPORT_cungqr F77_FUNC( cungqr , CUNGQR )
#define LAPACK_EXPORT_cungr2 F77_FUNC( cungr2 , CUNGR2 )
#define LAPACK_EXPORT_cungrq F77_FUNC( cungrq , CUNGRQ )
#define LAPACK_EXPORT_cungtr F77_FUNC( cungtr , CUNGTR )
#define LAPACK_EXPORT_cungtsqr F77_FUNC( cungtsqr , CUNGTSQR )
#define LAPACK_EXPORT_cunhr_col F77_FUNC( cunhr_col , CUNHR_COL )
#define LAPACK_EXPORT_cunm22 F77_FUNC( cunm22 , CUNM22 )
#define LAPACK_EXPORT_cunm2l F77_FUNC( cunm2l , CUNM2L )
#define LAPACK_EXPORT_cunm2r F77_FUNC( cunm2r , CUNM2R )
#define LAPACK_EXPORT_cunmbr F77_FUNC( cunmbr , CUNMBR )
#define LAPACK_EXPORT_cunmhr F77_FUNC( cunmhr , CUNMHR )
#define LAPACK_EXPORT_cunml2 F77_FUNC( cunml2 , CUNML2 )
#define LAPACK_EXPORT_cunmlq F77_FUNC( cunmlq , CUNMLQ )
#define LAPACK_EXPORT_cunmql F77_FUNC( cunmql , CUNMQL )
#define LAPACK_EXPORT_cunmqr F77_FUNC( cunmqr , CUNMQR )
#define LAPACK_EXPORT_cunmr2 F77_FUNC( cunmr2 , CUNMR2 )
#define LAPACK_EXPORT_cunmr3 F77_FUNC( cunmr3 , CUNMR3 )
#define LAPACK_EXPORT_cunmrq F77_FUNC( cunmrq , CUNMRQ )
#define LAPACK_EXPORT_cunmrz F77_FUNC( cunmrz , CUNMRZ )
#define LAPACK_EXPORT_cunmtr F77_FUNC( cunmtr , CUNMTR )
#define LAPACK_EXPORT_cupgtr F77_FUNC( cupgtr , CUPGTR )
#define LAPACK_EXPORT_cupmtr F77_FUNC( cupmtr , CUPMTR )
#define LAPACK_EXPORT_dbbcsd F77_FUNC( dbbcsd , DBBCSD )
#define LAPACK_EXPORT_dbdsdc F77_FUNC( dbdsdc , DBDSDC )
#define LAPACK_EXPORT_dbdsqr F77_FUNC( dbdsqr , DBDSQR )
#define LAPACK_EXPORT_dbdsvdx F77_FUNC( dbdsvdx , DBDSVDX )
#define LAPACK_EXPORT_dcombssq F77_FUNC( dcombssq , DCOMBSSQ )
#define LAPACK_EXPORT_ddisna F77_FUNC( ddisna , DDISNA )
#define LAPACK_EXPORT_dgbbrd F77_FUNC( dgbbrd , DGBBRD )
#define LAPACK_EXPORT_dgbcon F77_FUNC( dgbcon , DGBCON )
#define LAPACK_EXPORT_dgbequb F77_FUNC( dgbequb , DGBEQUB )
#define LAPACK_EXPORT_dgbequ F77_FUNC( dgbequ , DGBEQU )
#define LAPACK_EXPORT_dgbrfs F77_FUNC( dgbrfs , DGBRFS )
#define LAPACK_EXPORT_dgbrfsx F77_FUNC( dgbrfsx , DGBRFSX )
#define LAPACK_EXPORT_dgbsv F77_FUNC( dgbsv , DGBSV )
#define LAPACK_EXPORT_dgbsvx F77_FUNC( dgbsvx , DGBSVX )
#define LAPACK_EXPORT_dgbsvxx F77_FUNC( dgbsvxx , DGBSVXX )
#define LAPACK_EXPORT_dgbtf2 F77_FUNC( dgbtf2 , DGBTF2 )
#define LAPACK_EXPORT_dgbtrf F77_FUNC( dgbtrf , DGBTRF )
#define LAPACK_EXPORT_dgbtrs F77_FUNC( dgbtrs , DGBTRS )
#define LAPACK_EXPORT_dgebak F77_FUNC( dgebak , DGEBAK )
#define LAPACK_EXPORT_dgebal F77_FUNC( dgebal , DGEBAL )
#define LAPACK_EXPORT_dgebd2 F77_FUNC( dgebd2 , DGEBD2 )
#define LAPACK_EXPORT_dgebrd F77_FUNC( dgebrd , DGEBRD )
#define LAPACK_EXPORT_dgecon F77_FUNC( dgecon , DGECON )
#define LAPACK_EXPORT_dgeequb F77_FUNC( dgeequb , DGEEQUB )
#define LAPACK_EXPORT_dgeequ F77_FUNC( dgeequ , DGEEQU )
#define LAPACK_EXPORT_dgees F77_FUNC( dgees , DGEES )
#define LAPACK_EXPORT_dgeesx F77_FUNC( dgeesx , DGEESX )
#define LAPACK_EXPORT_dgeev F77_FUNC( dgeev , DGEEV )
#define LAPACK_EXPORT_dgeevx F77_FUNC( dgeevx , DGEEVX )
#define LAPACK_EXPORT_dgegs F77_FUNC( dgegs , DGEGS )
#define LAPACK_EXPORT_dgegv F77_FUNC( dgegv , DGEGV )
#define LAPACK_EXPORT_dgehd2 F77_FUNC( dgehd2 , DGEHD2 )
#define LAPACK_EXPORT_dgehrd F77_FUNC( dgehrd , DGEHRD )
#define LAPACK_EXPORT_dgejsv F77_FUNC( dgejsv , DGEJSV )
#define LAPACK_EXPORT_dgelq2 F77_FUNC( dgelq2 , DGELQ2 )
#define LAPACK_EXPORT_dgelq F77_FUNC( dgelq , DGELQ )
#define LAPACK_EXPORT_dgelqf F77_FUNC( dgelqf , DGELQF )
#define LAPACK_EXPORT_dgelqt3 F77_FUNC( dgelqt3 , DGELQT3 )
#define LAPACK_EXPORT_dgelqt F77_FUNC( dgelqt , DGELQT )
#define LAPACK_EXPORT_dgelsd F77_FUNC( dgelsd , DGELSD )
#define LAPACK_EXPORT_dgels F77_FUNC( dgels , DGELS )
#define LAPACK_EXPORT_dgelss F77_FUNC( dgelss , DGELSS )
#define LAPACK_EXPORT_dgelsx F77_FUNC( dgelsx , DGELSX )
#define LAPACK_EXPORT_dgelsy F77_FUNC( dgelsy , DGELSY )
#define LAPACK_EXPORT_dgemlq F77_FUNC( dgemlq , DGEMLQ )
#define LAPACK_EXPORT_dgemlqt F77_FUNC( dgemlqt , DGEMLQT )
#define LAPACK_EXPORT_dgemqr F77_FUNC( dgemqr , DGEMQR )
#define LAPACK_EXPORT_dgemqrt F77_FUNC( dgemqrt , DGEMQRT )
#define LAPACK_EXPORT_dgeql2 F77_FUNC( dgeql2 , DGEQL2 )
#define LAPACK_EXPORT_dgeqlf F77_FUNC( dgeqlf , DGEQLF )
#define LAPACK_EXPORT_dgeqp3 F77_FUNC( dgeqp3 , DGEQP3 )
#define LAPACK_EXPORT_dgeqpf F77_FUNC( dgeqpf , DGEQPF )
#define LAPACK_EXPORT_dgeqr2 F77_FUNC( dgeqr2 , DGEQR2 )
#define LAPACK_EXPORT_dgeqr2p F77_FUNC( dgeqr2p , DGEQR2P )
#define LAPACK_EXPORT_dgeqr F77_FUNC( dgeqr , DGEQR )
#define LAPACK_EXPORT_dgeqrf F77_FUNC( dgeqrf , DGEQRF )
#define LAPACK_EXPORT_dgeqrfp F77_FUNC( dgeqrfp , DGEQRFP )
#define LAPACK_EXPORT_dgeqrt2 F77_FUNC( dgeqrt2 , DGEQRT2 )
#define LAPACK_EXPORT_dgeqrt3 F77_FUNC( dgeqrt3 , DGEQRT3 )
#define LAPACK_EXPORT_dgeqrt F77_FUNC( dgeqrt , DGEQRT )
#define LAPACK_EXPORT_dgerfs F77_FUNC( dgerfs , DGERFS )
#define LAPACK_EXPORT_dgerfsx F77_FUNC( dgerfsx , DGERFSX )
#define LAPACK_EXPORT_dgerq2 F77_FUNC( dgerq2 , DGERQ2 )
#define LAPACK_EXPORT_dgerqf F77_FUNC( dgerqf , DGERQF )
#define LAPACK_EXPORT_dgesc2 F77_FUNC( dgesc2 , DGESC2 )
#define LAPACK_EXPORT_dgesdd F77_FUNC( dgesdd , DGESDD )
#define LAPACK_EXPORT_dgesvd F77_FUNC( dgesvd , DGESVD )
#define LAPACK_EXPORT_dgesvdq F77_FUNC( dgesvdq , DGESVDQ )
#define LAPACK_EXPORT_dgesvdx F77_FUNC( dgesvdx , DGESVDX )
#define LAPACK_EXPORT_dgesv F77_FUNC( dgesv , DGESV )
#define LAPACK_EXPORT_dgesvj F77_FUNC( dgesvj , DGESVJ )
#define LAPACK_EXPORT_dgesvx F77_FUNC( dgesvx , DGESVX )
#define LAPACK_EXPORT_dgesvxx F77_FUNC( dgesvxx , DGESVXX )
#define LAPACK_EXPORT_dgetc2 F77_FUNC( dgetc2 , DGETC2 )
#define LAPACK_EXPORT_dgetf2 F77_FUNC( dgetf2 , DGETF2 )
#define LAPACK_EXPORT_dgetrf2 F77_FUNC( dgetrf2 , DGETRF2 )
#define LAPACK_EXPORT_dgetrf F77_FUNC( dgetrf , DGETRF )
#define LAPACK_EXPORT_dgetri F77_FUNC( dgetri , DGETRI )
#define LAPACK_EXPORT_dgetrs F77_FUNC( dgetrs , DGETRS )
#define LAPACK_EXPORT_dgetsls F77_FUNC( dgetsls , DGETSLS )
#define LAPACK_EXPORT_dggbak F77_FUNC( dggbak , DGGBAK )
#define LAPACK_EXPORT_dggbal F77_FUNC( dggbal , DGGBAL )
#define LAPACK_EXPORT_dgges3 F77_FUNC( dgges3 , DGGES3 )
#define LAPACK_EXPORT_dgges F77_FUNC( dgges , DGGES )
#define LAPACK_EXPORT_dggesx F77_FUNC( dggesx , DGGESX )
#define LAPACK_EXPORT_dggev3 F77_FUNC( dggev3 , DGGEV3 )
#define LAPACK_EXPORT_dggev F77_FUNC( dggev , DGGEV )
#define LAPACK_EXPORT_dggevx F77_FUNC( dggevx , DGGEVX )
#define LAPACK_EXPORT_dggglm F77_FUNC( dggglm , DGGGLM )
#define LAPACK_EXPORT_dgghd3 F77_FUNC( dgghd3 , DGGHD3 )
#define LAPACK_EXPORT_dgghrd F77_FUNC( dgghrd , DGGHRD )
#define LAPACK_EXPORT_dgglse F77_FUNC( dgglse , DGGLSE )
#define LAPACK_EXPORT_dggqrf F77_FUNC( dggqrf , DGGQRF )
#define LAPACK_EXPORT_dggrqf F77_FUNC( dggrqf , DGGRQF )
#define LAPACK_EXPORT_dggsvd3 F77_FUNC( dggsvd3 , DGGSVD3 )
#define LAPACK_EXPORT_dggsvd F77_FUNC( dggsvd , DGGSVD )
#define LAPACK_EXPORT_dggsvp3 F77_FUNC( dggsvp3 , DGGSVP3 )
#define LAPACK_EXPORT_dggsvp F77_FUNC( dggsvp , DGGSVP )
#define LAPACK_EXPORT_dgsvj0 F77_FUNC( dgsvj0 , DGSVJ0 )
#define LAPACK_EXPORT_dgsvj1 F77_FUNC( dgsvj1 , DGSVJ1 )
#define LAPACK_EXPORT_dgtcon F77_FUNC( dgtcon , DGTCON )
#define LAPACK_EXPORT_dgtrfs F77_FUNC( dgtrfs , DGTRFS )
#define LAPACK_EXPORT_dgtsv F77_FUNC( dgtsv , DGTSV )
#define LAPACK_EXPORT_dgtsvx F77_FUNC( dgtsvx , DGTSVX )
#define LAPACK_EXPORT_dgttrf F77_FUNC( dgttrf , DGTTRF )
#define LAPACK_EXPORT_dgttrs F77_FUNC( dgttrs , DGTTRS )
#define LAPACK_EXPORT_dgtts2 F77_FUNC( dgtts2 , DGTTS2 )
#define LAPACK_EXPORT_dhgeqz F77_FUNC( dhgeqz , DHGEQZ )
#define LAPACK_EXPORT_dhsein F77_FUNC( dhsein , DHSEIN )
#define LAPACK_EXPORT_dhseqr F77_FUNC( dhseqr , DHSEQR )
#define LAPACK_EXPORT_disnan F77_FUNC( disnan , DISNAN )
#define LAPACK_EXPORT_dlabad F77_FUNC( dlabad , DLABAD )
#define LAPACK_EXPORT_dlabrd F77_FUNC( dlabrd , DLABRD )
#define LAPACK_EXPORT_dlacn2 F77_FUNC( dlacn2 , DLACN2 )
#define LAPACK_EXPORT_dlacon F77_FUNC( dlacon , DLACON )
#define LAPACK_EXPORT_dlacpy F77_FUNC( dlacpy , DLACPY )
#define LAPACK_EXPORT_dladiv F77_FUNC( dladiv , DLADIV )
#define LAPACK_EXPORT_dlae2 F77_FUNC( dlae2 , DLAE2 )
#define LAPACK_EXPORT_dlaebz F77_FUNC( dlaebz , DLAEBZ )
#define LAPACK_EXPORT_dlaed0 F77_FUNC( dlaed0 , DLAED0 )
#define LAPACK_EXPORT_dlaed1 F77_FUNC( dlaed1 , DLAED1 )
#define LAPACK_EXPORT_dlaed2 F77_FUNC( dlaed2 , DLAED2 )
#define LAPACK_EXPORT_dlaed3 F77_FUNC( dlaed3 , DLAED3 )
#define LAPACK_EXPORT_dlaed4 F77_FUNC( dlaed4 , DLAED4 )
#define LAPACK_EXPORT_dlaed5 F77_FUNC( dlaed5 , DLAED5 )
#define LAPACK_EXPORT_dlaed6 F77_FUNC( dlaed6 , DLAED6 )
#define LAPACK_EXPORT_dlaed7 F77_FUNC( dlaed7 , DLAED7 )
#define LAPACK_EXPORT_dlaed8 F77_FUNC( dlaed8 , DLAED8 )
#define LAPACK_EXPORT_dlaed9 F77_FUNC( dlaed9 , DLAED9 )
#define LAPACK_EXPORT_dlaeda F77_FUNC( dlaeda , DLAEDA )
#define LAPACK_EXPORT_dlaein F77_FUNC( dlaein , DLAEIN )
#define LAPACK_EXPORT_dlaev2 F77_FUNC( dlaev2 , DLAEV2 )
#define LAPACK_EXPORT_dlaexc F77_FUNC( dlaexc , DLAEXC )
#define LAPACK_EXPORT_dlag2 F77_FUNC( dlag2 , DLAG2 )
#define LAPACK_EXPORT_dlag2s F77_FUNC( dlag2s , DLAG2S )
#define LAPACK_EXPORT_dla_gbamv F77_FUNC( dla_gbamv , DLA_GBAMV )
#define LAPACK_EXPORT_dla_gbrcond F77_FUNC( dla_gbrcond , DLA_GBRCOND )
#define LAPACK_EXPORT_dla_gbrfsx_extended F77_FUNC( dla_gbrfsx_extended , DLA_GBRFSX_EXTENDED )
#define LAPACK_EXPORT_dla_gbrpvgrw F77_FUNC( dla_gbrpvgrw , DLA_GBRPVGRW )
#define LAPACK_EXPORT_dla_geamv F77_FUNC( dla_geamv , DLA_GEAMV )
#define LAPACK_EXPORT_dla_gercond F77_FUNC( dla_gercond , DLA_GERCOND )
#define LAPACK_EXPORT_dla_gerfsx_extended F77_FUNC( dla_gerfsx_extended , DLA_GERFSX_EXTENDED )
#define LAPACK_EXPORT_dla_gerpvgrw F77_FUNC( dla_gerpvgrw , DLA_GERPVGRW )
#define LAPACK_EXPORT_dlags2 F77_FUNC( dlags2 , DLAGS2 )
#define LAPACK_EXPORT_dlagtf F77_FUNC( dlagtf , DLAGTF )
#define LAPACK_EXPORT_dlagtm F77_FUNC( dlagtm , DLAGTM )
#define LAPACK_EXPORT_dlagts F77_FUNC( dlagts , DLAGTS )
#define LAPACK_EXPORT_dlagv2 F77_FUNC( dlagv2 , DLAGV2 )
#define LAPACK_EXPORT_dlahqr F77_FUNC( dlahqr , DLAHQR )
#define LAPACK_EXPORT_dlahr2 F77_FUNC( dlahr2 , DLAHR2 )
#define LAPACK_EXPORT_dlahrd F77_FUNC( dlahrd , DLAHRD )
#define LAPACK_EXPORT_dlaic1 F77_FUNC( dlaic1 , DLAIC1 )
#define LAPACK_EXPORT_dlaisnan F77_FUNC( dlaisnan , DLAISNAN )
#define LAPACK_EXPORT_dla_lin_berr F77_FUNC( dla_lin_berr , DLA_LIN_BERR )
#define LAPACK_EXPORT_dlaln2 F77_FUNC( dlaln2 , DLALN2 )
#define LAPACK_EXPORT_dlals0 F77_FUNC( dlals0 , DLALS0 )
#define LAPACK_EXPORT_dlalsa F77_FUNC( dlalsa , DLALSA )
#define LAPACK_EXPORT_dlalsd F77_FUNC( dlalsd , DLALSD )
#define LAPACK_EXPORT_dlamrg F77_FUNC( dlamrg , DLAMRG )
#define LAPACK_EXPORT_dlamswlq F77_FUNC( dlamswlq , DLAMSWLQ )
#define LAPACK_EXPORT_dlamtsqr F77_FUNC( dlamtsqr , DLAMTSQR )
#define LAPACK_EXPORT_dlaneg F77_FUNC( dlaneg , DLANEG )
#define LAPACK_EXPORT_dlangb F77_FUNC( dlangb , DLANGB )
#define LAPACK_EXPORT_dlange F77_FUNC( dlange , DLANGE )
#define LAPACK_EXPORT_dlangt F77_FUNC( dlangt , DLANGT )
#define LAPACK_EXPORT_dlanhs F77_FUNC( dlanhs , DLANHS )
#define LAPACK_EXPORT_dlansb F77_FUNC( dlansb , DLANSB )
#define LAPACK_EXPORT_dlansf F77_FUNC( dlansf , DLANSF )
#define LAPACK_EXPORT_dlansp F77_FUNC( dlansp , DLANSP )
#define LAPACK_EXPORT_dlanst F77_FUNC( dlanst , DLANST )
#define LAPACK_EXPORT_dlansy F77_FUNC( dlansy , DLANSY )
#define LAPACK_EXPORT_dlantb F77_FUNC( dlantb , DLANTB )
#define LAPACK_EXPORT_dlantp F77_FUNC( dlantp , DLANTP )
#define LAPACK_EXPORT_dlantr F77_FUNC( dlantr , DLANTR )
#define LAPACK_EXPORT_dlanv2 F77_FUNC( dlanv2 , DLANV2 )
#define LAPACK_EXPORT_dlaorhr_col_getrfnp2 F77_FUNC( dlaorhr_col_getrfnp2 , DLAORHR_COL_GETRFNP2 )
#define LAPACK_EXPORT_dlaorhr_col_getrfnp F77_FUNC( dlaorhr_col_getrfnp , DLAORHR_COL_GETRFNP )
#define LAPACK_EXPORT_dlapll F77_FUNC( dlapll , DLAPLL )
#define LAPACK_EXPORT_dlapmr F77_FUNC( dlapmr , DLAPMR )
#define LAPACK_EXPORT_dlapmt F77_FUNC( dlapmt , DLAPMT )
#define LAPACK_EXPORT_dla_porcond F77_FUNC( dla_porcond , DLA_PORCOND )
#define LAPACK_EXPORT_dla_porfsx_extended F77_FUNC( dla_porfsx_extended , DLA_PORFSX_EXTENDED )
#define LAPACK_EXPORT_dla_porpvgrw F77_FUNC( dla_porpvgrw , DLA_PORPVGRW )
#define LAPACK_EXPORT_dlapy2 F77_FUNC( dlapy2 , DLAPY2 )
#define LAPACK_EXPORT_dlapy3 F77_FUNC( dlapy3 , DLAPY3 )
#define LAPACK_EXPORT_dlaqgb F77_FUNC( dlaqgb , DLAQGB )
#define LAPACK_EXPORT_dlaqge F77_FUNC( dlaqge , DLAQGE )
#define LAPACK_EXPORT_dlaqp2 F77_FUNC( dlaqp2 , DLAQP2 )
#define LAPACK_EXPORT_dlaqps F77_FUNC( dlaqps , DLAQPS )
#define LAPACK_EXPORT_dlaqr0 F77_FUNC( dlaqr0 , DLAQR0 )
#define LAPACK_EXPORT_dlaqr1 F77_FUNC( dlaqr1 , DLAQR1 )
#define LAPACK_EXPORT_dlaqr2 F77_FUNC( dlaqr2 , DLAQR2 )
#define LAPACK_EXPORT_dlaqr3 F77_FUNC( dlaqr3 , DLAQR3 )
#define LAPACK_EXPORT_dlaqr4 F77_FUNC( dlaqr4 , DLAQR4 )
#define LAPACK_EXPORT_dlaqr5 F77_FUNC( dlaqr5 , DLAQR5 )
#define LAPACK_EXPORT_dlaqsb F77_FUNC( dlaqsb , DLAQSB )
#define LAPACK_EXPORT_dlaqsp F77_FUNC( dlaqsp , DLAQSP )
#define LAPACK_EXPORT_dlaqsy F77_FUNC( dlaqsy , DLAQSY )
#define LAPACK_EXPORT_dlaqtr F77_FUNC( dlaqtr , DLAQTR )
#define LAPACK_EXPORT_dlar1v F77_FUNC( dlar1v , DLAR1V )
#define LAPACK_EXPORT_dlar2v F77_FUNC( dlar2v , DLAR2V )
#define LAPACK_EXPORT_dlarfb F77_FUNC( dlarfb , DLARFB )
#define LAPACK_EXPORT_dlarf F77_FUNC( dlarf , DLARF )
#define LAPACK_EXPORT_dlarfg F77_FUNC( dlarfg , DLARFG )
#define LAPACK_EXPORT_dlarfgp F77_FUNC( dlarfgp , DLARFGP )
#define LAPACK_EXPORT_dlarft F77_FUNC( dlarft , DLARFT )
#define LAPACK_EXPORT_dlarfx F77_FUNC( dlarfx , DLARFX )
#define LAPACK_EXPORT_dlarfy F77_FUNC( dlarfy , DLARFY )
#define LAPACK_EXPORT_dlargv F77_FUNC( dlargv , DLARGV )
#define LAPACK_EXPORT_dlarnv F77_FUNC( dlarnv , DLARNV )
#define LAPACK_EXPORT_dlarra F77_FUNC( dlarra , DLARRA )
#define LAPACK_EXPORT_dlarrb F77_FUNC( dlarrb , DLARRB )
#define LAPACK_EXPORT_dlarrc F77_FUNC( dlarrc , DLARRC )
#define LAPACK_EXPORT_dlarrd F77_FUNC( dlarrd , DLARRD )
#define LAPACK_EXPORT_dlarre F77_FUNC( dlarre , DLARRE )
#define LAPACK_EXPORT_dlarrf F77_FUNC( dlarrf , DLARRF )
#define LAPACK_EXPORT_dlarrj F77_FUNC( dlarrj , DLARRJ )
#define LAPACK_EXPORT_dlarrk F77_FUNC( dlarrk , DLARRK )
#define LAPACK_EXPORT_dlarrr F77_FUNC( dlarrr , DLARRR )
#define LAPACK_EXPORT_dlarrv F77_FUNC( dlarrv , DLARRV )
#define LAPACK_EXPORT_dlarscl2 F77_FUNC( dlarscl2 , DLARSCL2 )
#define LAPACK_EXPORT_dlartg F77_FUNC( dlartg , DLARTG )
#define LAPACK_EXPORT_dlartgp F77_FUNC( dlartgp , DLARTGP )
#define LAPACK_EXPORT_dlartgs F77_FUNC( dlartgs , DLARTGS )
#define LAPACK_EXPORT_dlartv F77_FUNC( dlartv , DLARTV )
#define LAPACK_EXPORT_dlaruv F77_FUNC( dlaruv , DLARUV )
#define LAPACK_EXPORT_dlarzb F77_FUNC( dlarzb , DLARZB )
#define LAPACK_EXPORT_dlarz F77_FUNC( dlarz , DLARZ )
#define LAPACK_EXPORT_dlarzt F77_FUNC( dlarzt , DLARZT )
#define LAPACK_EXPORT_dlas2 F77_FUNC( dlas2 , DLAS2 )
#define LAPACK_EXPORT_dlascl2 F77_FUNC( dlascl2 , DLASCL2 )
#define LAPACK_EXPORT_dlascl F77_FUNC( dlascl , DLASCL )
#define LAPACK_EXPORT_dlasd0 F77_FUNC( dlasd0 , DLASD0 )
#define LAPACK_EXPORT_dlasd1 F77_FUNC( dlasd1 , DLASD1 )
#define LAPACK_EXPORT_dlasd2 F77_FUNC( dlasd2 , DLASD2 )
#define LAPACK_EXPORT_dlasd3 F77_FUNC( dlasd3 , DLASD3 )
#define LAPACK_EXPORT_dlasd4 F77_FUNC( dlasd4 , DLASD4 )
#define LAPACK_EXPORT_dlasd5 F77_FUNC( dlasd5 , DLASD5 )
#define LAPACK_EXPORT_dlasd6 F77_FUNC( dlasd6 , DLASD6 )
#define LAPACK_EXPORT_dlasd7 F77_FUNC( dlasd7 , DLASD7 )
#define LAPACK_EXPORT_dlasd8 F77_FUNC( dlasd8 , DLASD8 )
#define LAPACK_EXPORT_dlasda F77_FUNC( dlasda , DLASDA )
#define LAPACK_EXPORT_dlasdq F77_FUNC( dlasdq , DLASDQ )
#define LAPACK_EXPORT_dlasdt F77_FUNC( dlasdt , DLASDT )
#define LAPACK_EXPORT_dlaset F77_FUNC( dlaset , DLASET )
#define LAPACK_EXPORT_dlasq1 F77_FUNC( dlasq1 , DLASQ1 )
#define LAPACK_EXPORT_dlasq2 F77_FUNC( dlasq2 , DLASQ2 )
#define LAPACK_EXPORT_dlasq3 F77_FUNC( dlasq3 , DLASQ3 )
#define LAPACK_EXPORT_dlasq4 F77_FUNC( dlasq4 , DLASQ4 )
#define LAPACK_EXPORT_dlasq5 F77_FUNC( dlasq5 , DLASQ5 )
#define LAPACK_EXPORT_dlasq6 F77_FUNC( dlasq6 , DLASQ6 )
#define LAPACK_EXPORT_dlasr F77_FUNC( dlasr , DLASR )
#define LAPACK_EXPORT_dlasrt F77_FUNC( dlasrt , DLASRT )
#define LAPACK_EXPORT_dlassq F77_FUNC( dlassq , DLASSQ )
#define LAPACK_EXPORT_dlasv2 F77_FUNC( dlasv2 , DLASV2 )
#define LAPACK_EXPORT_dlaswlq F77_FUNC( dlaswlq , DLASWLQ )
#define LAPACK_EXPORT_dlaswp F77_FUNC( dlaswp , DLASWP )
#define LAPACK_EXPORT_dlasy2 F77_FUNC( dlasy2 , DLASY2 )
#define LAPACK_EXPORT_dla_syamv F77_FUNC( dla_syamv , DLA_SYAMV )
#define LAPACK_EXPORT_dlasyf_aa F77_FUNC( dlasyf_aa , DLASYF_AA )
#define LAPACK_EXPORT_dlasyf F77_FUNC( dlasyf , DLASYF )
#define LAPACK_EXPORT_dlasyf_rk F77_FUNC( dlasyf_rk , DLASYF_RK )
#define LAPACK_EXPORT_dlasyf_rook F77_FUNC( dlasyf_rook , DLASYF_ROOK )
#define LAPACK_EXPORT_dla_syrcond F77_FUNC( dla_syrcond , DLA_SYRCOND )
#define LAPACK_EXPORT_dla_syrfsx_extended F77_FUNC( dla_syrfsx_extended , DLA_SYRFSX_EXTENDED )
#define LAPACK_EXPORT_dla_syrpvgrw F77_FUNC( dla_syrpvgrw , DLA_SYRPVGRW )
#define LAPACK_EXPORT_dlat2s F77_FUNC( dlat2s , DLAT2S )
#define LAPACK_EXPORT_dlatbs F77_FUNC( dlatbs , DLATBS )
#define LAPACK_EXPORT_dlatdf F77_FUNC( dlatdf , DLATDF )
#define LAPACK_EXPORT_dlatps F77_FUNC( dlatps , DLATPS )
#define LAPACK_EXPORT_dlatrd F77_FUNC( dlatrd , DLATRD )
#define LAPACK_EXPORT_dlatrs F77_FUNC( dlatrs , DLATRS )
#define LAPACK_EXPORT_dlatrz F77_FUNC( dlatrz , DLATRZ )
#define LAPACK_EXPORT_dlatsqr F77_FUNC( dlatsqr , DLATSQR )
#define LAPACK_EXPORT_dlatzm F77_FUNC( dlatzm , DLATZM )
#define LAPACK_EXPORT_dlauu2 F77_FUNC( dlauu2 , DLAUU2 )
#define LAPACK_EXPORT_dlauum F77_FUNC( dlauum , DLAUUM )
#define LAPACK_EXPORT_dla_wwaddw F77_FUNC( dla_wwaddw , DLA_WWADDW )
#define LAPACK_EXPORT_dopmtr F77_FUNC( dopmtr , DOPMTR )
#define LAPACK_EXPORT_dorbdb1 F77_FUNC( dorbdb1 , DORBDB1 )
#define LAPACK_EXPORT_dorbdb2 F77_FUNC( dorbdb2 , DORBDB2 )
#define LAPACK_EXPORT_dorbdb3 F77_FUNC( dorbdb3 , DORBDB3 )
#define LAPACK_EXPORT_dorbdb4 F77_FUNC( dorbdb4 , DORBDB4 )
#define LAPACK_EXPORT_dorbdb5 F77_FUNC( dorbdb5 , DORBDB5 )
#define LAPACK_EXPORT_dorbdb6 F77_FUNC( dorbdb6 , DORBDB6 )
#define LAPACK_EXPORT_dorbdb F77_FUNC( dorbdb , DORBDB )
#define LAPACK_EXPORT_dorg2l F77_FUNC( dorg2l , DORG2L )
#define LAPACK_EXPORT_dorg2r F77_FUNC( dorg2r , DORG2R )
#define LAPACK_EXPORT_dorgbr F77_FUNC( dorgbr , DORGBR )
#define LAPACK_EXPORT_dorgl2 F77_FUNC( dorgl2 , DORGL2 )
#define LAPACK_EXPORT_dorglq F77_FUNC( dorglq , DORGLQ )
#define LAPACK_EXPORT_dorgql F77_FUNC( dorgql , DORGQL )
#define LAPACK_EXPORT_dorgqr F77_FUNC( dorgqr , DORGQR )
#define LAPACK_EXPORT_dorgr2 F77_FUNC( dorgr2 , DORGR2 )
#define LAPACK_EXPORT_dorgrq F77_FUNC( dorgrq , DORGRQ )
#define LAPACK_EXPORT_dorgtr F77_FUNC( dorgtr , DORGTR )
#define LAPACK_EXPORT_dorgtsqr F77_FUNC( dorgtsqr , DORGTSQR )
#define LAPACK_EXPORT_dorhr_col F77_FUNC( dorhr_col , DORHR_COL )
#define LAPACK_EXPORT_dorm22 F77_FUNC( dorm22 , DORM22 )
#define LAPACK_EXPORT_dorm2l F77_FUNC( dorm2l , DORM2L )
#define LAPACK_EXPORT_dorm2r F77_FUNC( dorm2r , DORM2R )
#define LAPACK_EXPORT_dormbr F77_FUNC( dormbr , DORMBR )
#define LAPACK_EXPORT_dorml2 F77_FUNC( dorml2 , DORML2 )
#define LAPACK_EXPORT_dormlq F77_FUNC( dormlq , DORMLQ )
#define LAPACK_EXPORT_dormql F77_FUNC( dormql , DORMQL )
#define LAPACK_EXPORT_dormqr F77_FUNC( dormqr , DORMQR )
#define LAPACK_EXPORT_dormr2 F77_FUNC( dormr2 , DORMR2 )
#define LAPACK_EXPORT_dormr3 F77_FUNC( dormr3 , DORMR3 )
#define LAPACK_EXPORT_dormrq F77_FUNC( dormrq , DORMRQ )
#define LAPACK_EXPORT_dormrz F77_FUNC( dormrz , DORMRZ )
#define LAPACK_EXPORT_dormtr F77_FUNC( dormtr , DORMTR )
#define LAPACK_EXPORT_dpbcon F77_FUNC( dpbcon , DPBCON )
#define LAPACK_EXPORT_dpbequ F77_FUNC( dpbequ , DPBEQU )
#define LAPACK_EXPORT_dpbrfs F77_FUNC( dpbrfs , DPBRFS )
#define LAPACK_EXPORT_dpbstf F77_FUNC( dpbstf , DPBSTF )
#define LAPACK_EXPORT_dpbsv F77_FUNC( dpbsv , DPBSV )
#define LAPACK_EXPORT_dpbsvx F77_FUNC( dpbsvx , DPBSVX )
#define LAPACK_EXPORT_dpbtf2 F77_FUNC( dpbtf2 , DPBTF2 )
#define LAPACK_EXPORT_dpbtrf F77_FUNC( dpbtrf , DPBTRF )
#define LAPACK_EXPORT_dpbtrs F77_FUNC( dpbtrs , DPBTRS )
#define LAPACK_EXPORT_dpftrf F77_FUNC( dpftrf , DPFTRF )
#define LAPACK_EXPORT_dpftri F77_FUNC( dpftri , DPFTRI )
#define LAPACK_EXPORT_dpftrs F77_FUNC( dpftrs , DPFTRS )
#define LAPACK_EXPORT_dpocon F77_FUNC( dpocon , DPOCON )
#define LAPACK_EXPORT_dpoequb F77_FUNC( dpoequb , DPOEQUB )
#define LAPACK_EXPORT_dpoequ F77_FUNC( dpoequ , DPOEQU )
#define LAPACK_EXPORT_dporfs F77_FUNC( dporfs , DPORFS )
#define LAPACK_EXPORT_dporfsx F77_FUNC( dporfsx , DPORFSX )
#define LAPACK_EXPORT_dposv F77_FUNC( dposv , DPOSV )
#define LAPACK_EXPORT_dposvx F77_FUNC( dposvx , DPOSVX )
#define LAPACK_EXPORT_dposvxx F77_FUNC( dposvxx , DPOSVXX )
#define LAPACK_EXPORT_dpotf2 F77_FUNC( dpotf2 , DPOTF2 )
#define LAPACK_EXPORT_dpotrf2 F77_FUNC( dpotrf2 , DPOTRF2 )
#define LAPACK_EXPORT_dpotrf F77_FUNC( dpotrf , DPOTRF )
#define LAPACK_EXPORT_dpotri F77_FUNC( dpotri , DPOTRI )
#define LAPACK_EXPORT_dpotrs F77_FUNC( dpotrs , DPOTRS )
#define LAPACK_EXPORT_dppcon F77_FUNC( dppcon , DPPCON )
#define LAPACK_EXPORT_dppequ F77_FUNC( dppequ , DPPEQU )
#define LAPACK_EXPORT_dpprfs F77_FUNC( dpprfs , DPPRFS )
#define LAPACK_EXPORT_dppsv F77_FUNC( dppsv , DPPSV )
#define LAPACK_EXPORT_dppsvx F77_FUNC( dppsvx , DPPSVX )
#define LAPACK_EXPORT_dpptrf F77_FUNC( dpptrf , DPPTRF )
#define LAPACK_EXPORT_dpptri F77_FUNC( dpptri , DPPTRI )
#define LAPACK_EXPORT_dpptrs F77_FUNC( dpptrs , DPPTRS )
#define LAPACK_EXPORT_dpstf2 F77_FUNC( dpstf2 , DPSTF2 )
#define LAPACK_EXPORT_dpstrf F77_FUNC( dpstrf , DPSTRF )
#define LAPACK_EXPORT_dptcon F77_FUNC( dptcon , DPTCON )
#define LAPACK_EXPORT_dpteqr F77_FUNC( dpteqr , DPTEQR )
#define LAPACK_EXPORT_dptrfs F77_FUNC( dptrfs , DPTRFS )
#define LAPACK_EXPORT_dptsv F77_FUNC( dptsv , DPTSV )
#define LAPACK_EXPORT_dptsvx F77_FUNC( dptsvx , DPTSVX )
#define LAPACK_EXPORT_dpttrf F77_FUNC( dpttrf , DPTTRF )
#define LAPACK_EXPORT_dpttrs F77_FUNC( dpttrs , DPTTRS )
#define LAPACK_EXPORT_dptts2 F77_FUNC( dptts2 , DPTTS2 )
#define LAPACK_EXPORT_drscl F77_FUNC( drscl , DRSCL )
#define LAPACK_EXPORT_dsb2st_kernels F77_FUNC( dsb2st_kernels , DSB2ST_KERNELS )
#define LAPACK_EXPORT_dsbev_2stage F77_FUNC( dsbev_2stage , DSBEV_2STAGE )
#define LAPACK_EXPORT_dsbevd_2stage F77_FUNC( dsbevd_2stage , DSBEVD_2STAGE )
#define LAPACK_EXPORT_dsbevd F77_FUNC( dsbevd , DSBEVD )
#define LAPACK_EXPORT_dsbev F77_FUNC( dsbev , DSBEV )
#define LAPACK_EXPORT_dsbevx_2stage F77_FUNC( dsbevx_2stage , DSBEVX_2STAGE )
#define LAPACK_EXPORT_dsbevx F77_FUNC( dsbevx , DSBEVX )
#define LAPACK_EXPORT_dsbgst F77_FUNC( dsbgst , DSBGST )
#define LAPACK_EXPORT_dsbgvd F77_FUNC( dsbgvd , DSBGVD )
#define LAPACK_EXPORT_dsbgv F77_FUNC( dsbgv , DSBGV )
#define LAPACK_EXPORT_dsbgvx F77_FUNC( dsbgvx , DSBGVX )
#define LAPACK_EXPORT_dsbtrd F77_FUNC( dsbtrd , DSBTRD )
#define LAPACK_EXPORT_dsfrk F77_FUNC( dsfrk , DSFRK )
#define LAPACK_EXPORT_dsgesv F77_FUNC( dsgesv , DSGESV )
#define LAPACK_EXPORT_dspcon F77_FUNC( dspcon , DSPCON )
#define LAPACK_EXPORT_dspevd F77_FUNC( dspevd , DSPEVD )
#define LAPACK_EXPORT_dspev F77_FUNC( dspev , DSPEV )
#define LAPACK_EXPORT_dspevx F77_FUNC( dspevx , DSPEVX )
#define LAPACK_EXPORT_dspgst F77_FUNC( dspgst , DSPGST )
#define LAPACK_EXPORT_dspgvd F77_FUNC( dspgvd , DSPGVD )
#define LAPACK_EXPORT_dspgv F77_FUNC( dspgv , DSPGV )
#define LAPACK_EXPORT_dspgvx F77_FUNC( dspgvx , DSPGVX )
#define LAPACK_EXPORT_dsposv F77_FUNC( dsposv , DSPOSV )
#define LAPACK_EXPORT_dsprfs F77_FUNC( dsprfs , DSPRFS )
#define LAPACK_EXPORT_dspsv F77_FUNC( dspsv , DSPSV )
#define LAPACK_EXPORT_dspsvx F77_FUNC( dspsvx , DSPSVX )
#define LAPACK_EXPORT_dsptrd F77_FUNC( dsptrd , DSPTRD )
#define LAPACK_EXPORT_dsptrf F77_FUNC( dsptrf , DSPTRF )
#define LAPACK_EXPORT_dsptri F77_FUNC( dsptri , DSPTRI )
#define LAPACK_EXPORT_dsptrs F77_FUNC( dsptrs , DSPTRS )
#define LAPACK_EXPORT_dstebz F77_FUNC( dstebz , DSTEBZ )
#define LAPACK_EXPORT_dstedc F77_FUNC( dstedc , DSTEDC )
#define LAPACK_EXPORT_dstegr F77_FUNC( dstegr , DSTEGR )
#define LAPACK_EXPORT_dstein F77_FUNC( dstein , DSTEIN )
#define LAPACK_EXPORT_dstemr F77_FUNC( dstemr , DSTEMR )
#define LAPACK_EXPORT_dsteqr F77_FUNC( dsteqr , DSTEQR )
#define LAPACK_EXPORT_dsterf F77_FUNC( dsterf , DSTERF )
#define LAPACK_EXPORT_dstevd F77_FUNC( dstevd , DSTEVD )
#define LAPACK_EXPORT_dstev F77_FUNC( dstev , DSTEV )
#define LAPACK_EXPORT_dstevr F77_FUNC( dstevr , DSTEVR )
#define LAPACK_EXPORT_dstevx F77_FUNC( dstevx , DSTEVX )
#define LAPACK_EXPORT_dsycon_3 F77_FUNC( dsycon_3 , DSYCON_3 )
#define LAPACK_EXPORT_dsycon F77_FUNC( dsycon , DSYCON )
#define LAPACK_EXPORT_dsycon_rook F77_FUNC( dsycon_rook , DSYCON_ROOK )
#define LAPACK_EXPORT_dsyconv F77_FUNC( dsyconv , DSYCONV )
#define LAPACK_EXPORT_dsyconvf F77_FUNC( dsyconvf , DSYCONVF )
#define LAPACK_EXPORT_dsyconvf_rook F77_FUNC( dsyconvf_rook , DSYCONVF_ROOK )
#define LAPACK_EXPORT_dsyequb F77_FUNC( dsyequb , DSYEQUB )
#define LAPACK_EXPORT_dsyev_2stage F77_FUNC( dsyev_2stage , DSYEV_2STAGE )
#define LAPACK_EXPORT_dsyevd_2stage F77_FUNC( dsyevd_2stage , DSYEVD_2STAGE )
#define LAPACK_EXPORT_dsyevd F77_FUNC( dsyevd , DSYEVD )
#define LAPACK_EXPORT_dsyev  F77_FUNC( dsyev  , DSYEV  )
#define LAPACK_EXPORT_dsyev F77_FUNC( dsyev , DSYEV )
#define LAPACK_EXPORT_dsyevr_2stage F77_FUNC( dsyevr_2stage , DSYEVR_2STAGE )
#define LAPACK_EXPORT_dsyevr F77_FUNC( dsyevr , DSYEVR )
#define LAPACK_EXPORT_dsyevx_2stage F77_FUNC( dsyevx_2stage , DSYEVX_2STAGE )
#define LAPACK_EXPORT_dsyevx F77_FUNC( dsyevx , DSYEVX )
#define LAPACK_EXPORT_dsygs2 F77_FUNC( dsygs2 , DSYGS2 )
#define LAPACK_EXPORT_dsygst F77_FUNC( dsygst , DSYGST )
#define LAPACK_EXPORT_dsygv_2stage F77_FUNC( dsygv_2stage , DSYGV_2STAGE )
#define LAPACK_EXPORT_dsygvd F77_FUNC( dsygvd , DSYGVD )
#define LAPACK_EXPORT_dsygv F77_FUNC( dsygv , DSYGV )
#define LAPACK_EXPORT_dsygvx F77_FUNC( dsygvx , DSYGVX )
#define LAPACK_EXPORT_dsyrfs F77_FUNC( dsyrfs , DSYRFS )
#define LAPACK_EXPORT_dsyrfsx F77_FUNC( dsyrfsx , DSYRFSX )
#define LAPACK_EXPORT_dsysv_aa_2stage F77_FUNC( dsysv_aa_2stage , DSYSV_AA_2STAGE )
#define LAPACK_EXPORT_dsysv_aa F77_FUNC( dsysv_aa , DSYSV_AA )
#define LAPACK_EXPORT_dsysv F77_FUNC( dsysv , DSYSV )
#define LAPACK_EXPORT_dsysv_rk F77_FUNC( dsysv_rk , DSYSV_RK )
#define LAPACK_EXPORT_dsysv_rook F77_FUNC( dsysv_rook , DSYSV_ROOK )
#define LAPACK_EXPORT_dsysvx F77_FUNC( dsysvx , DSYSVX )
#define LAPACK_EXPORT_dsysvxx F77_FUNC( dsysvxx , DSYSVXX )
#define LAPACK_EXPORT_dsyswapr F77_FUNC( dsyswapr , DSYSWAPR )
#define LAPACK_EXPORT_dsytd2 F77_FUNC( dsytd2 , DSYTD2 )
#define LAPACK_EXPORT_dsytf2 F77_FUNC( dsytf2 , DSYTF2 )
#define LAPACK_EXPORT_dsytf2_rk F77_FUNC( dsytf2_rk , DSYTF2_RK )
#define LAPACK_EXPORT_dsytf2_rook F77_FUNC( dsytf2_rook , DSYTF2_ROOK )
#define LAPACK_EXPORT_dsytrd_2stage F77_FUNC( dsytrd_2stage , DSYTRD_2STAGE )
#define LAPACK_EXPORT_dsytrd F77_FUNC( dsytrd , DSYTRD )
#define LAPACK_EXPORT_dsytrd_sb2st F77_FUNC( dsytrd_sb2st , DSYTRD_SB2ST )
#define LAPACK_EXPORT_dsytrd_sy2sb F77_FUNC( dsytrd_sy2sb , DSYTRD_SY2SB )
#define LAPACK_EXPORT_dsytrf_aa_2stage F77_FUNC( dsytrf_aa_2stage , DSYTRF_AA_2STAGE )
#define LAPACK_EXPORT_dsytrf_aa F77_FUNC( dsytrf_aa , DSYTRF_AA )
#define LAPACK_EXPORT_dsytrf F77_FUNC( dsytrf , DSYTRF )
#define LAPACK_EXPORT_dsytrf_rk F77_FUNC( dsytrf_rk , DSYTRF_RK )
#define LAPACK_EXPORT_dsytrf_rook F77_FUNC( dsytrf_rook , DSYTRF_ROOK )
#define LAPACK_EXPORT_dsytri2 F77_FUNC( dsytri2 , DSYTRI2 )
#define LAPACK_EXPORT_dsytri2x F77_FUNC( dsytri2x , DSYTRI2X )
#define LAPACK_EXPORT_dsytri_3 F77_FUNC( dsytri_3 , DSYTRI_3 )
#define LAPACK_EXPORT_dsytri_3x F77_FUNC( dsytri_3x , DSYTRI_3X )
#define LAPACK_EXPORT_dsytri F77_FUNC( dsytri , DSYTRI )
#define LAPACK_EXPORT_dsytri_rook F77_FUNC( dsytri_rook , DSYTRI_ROOK )
#define LAPACK_EXPORT_dsytrs2 F77_FUNC( dsytrs2 , DSYTRS2 )
#define LAPACK_EXPORT_dsytrs_3 F77_FUNC( dsytrs_3 , DSYTRS_3 )
#define LAPACK_EXPORT_dsytrs_aa_2stage F77_FUNC( dsytrs_aa_2stage , DSYTRS_AA_2STAGE )
#define LAPACK_EXPORT_dsytrs_aa F77_FUNC( dsytrs_aa , DSYTRS_AA )
#define LAPACK_EXPORT_dsytrs F77_FUNC( dsytrs , DSYTRS )
#define LAPACK_EXPORT_dsytrs_rook F77_FUNC( dsytrs_rook , DSYTRS_ROOK )
#define LAPACK_EXPORT_dtbcon F77_FUNC( dtbcon , DTBCON )
#define LAPACK_EXPORT_dtbrfs F77_FUNC( dtbrfs , DTBRFS )
#define LAPACK_EXPORT_dtbtrs F77_FUNC( dtbtrs , DTBTRS )
#define LAPACK_EXPORT_dtfsm F77_FUNC( dtfsm , DTFSM )
#define LAPACK_EXPORT_dtftri F77_FUNC( dtftri , DTFTRI )
#define LAPACK_EXPORT_dtfttp F77_FUNC( dtfttp , DTFTTP )
#define LAPACK_EXPORT_dtfttr F77_FUNC( dtfttr , DTFTTR )
#define LAPACK_EXPORT_dtgevc F77_FUNC( dtgevc , DTGEVC )
#define LAPACK_EXPORT_dtgex2 F77_FUNC( dtgex2 , DTGEX2 )
#define LAPACK_EXPORT_dtgexc F77_FUNC( dtgexc , DTGEXC )
#define LAPACK_EXPORT_dtgsen F77_FUNC( dtgsen , DTGSEN )
#define LAPACK_EXPORT_dtgsja F77_FUNC( dtgsja , DTGSJA )
#define LAPACK_EXPORT_dtgsna F77_FUNC( dtgsna , DTGSNA )
#define LAPACK_EXPORT_dtgsy2 F77_FUNC( dtgsy2 , DTGSY2 )
#define LAPACK_EXPORT_dtgsyl F77_FUNC( dtgsyl , DTGSYL )
#define LAPACK_EXPORT_dtpcon F77_FUNC( dtpcon , DTPCON )
#define LAPACK_EXPORT_dtplqt2 F77_FUNC( dtplqt2 , DTPLQT2 )
#define LAPACK_EXPORT_dtplqt F77_FUNC( dtplqt , DTPLQT )
#define LAPACK_EXPORT_dtpmlqt F77_FUNC( dtpmlqt , DTPMLQT )
#define LAPACK_EXPORT_dtpmqrt F77_FUNC( dtpmqrt , DTPMQRT )
#define LAPACK_EXPORT_dtpqrt2 F77_FUNC( dtpqrt2 , DTPQRT2 )
#define LAPACK_EXPORT_dtpqrt F77_FUNC( dtpqrt , DTPQRT )
#define LAPACK_EXPORT_dtprfb F77_FUNC( dtprfb , DTPRFB )
#define LAPACK_EXPORT_dtprfs F77_FUNC( dtprfs , DTPRFS )
#define LAPACK_EXPORT_dtptri F77_FUNC( dtptri , DTPTRI )
#define LAPACK_EXPORT_dtptrs F77_FUNC( dtptrs , DTPTRS )
#define LAPACK_EXPORT_dtpttf F77_FUNC( dtpttf , DTPTTF )
#define LAPACK_EXPORT_dtpttr F77_FUNC( dtpttr , DTPTTR )
#define LAPACK_EXPORT_dtrcon F77_FUNC( dtrcon , DTRCON )
#define LAPACK_EXPORT_dtrevc3 F77_FUNC( dtrevc3 , DTREVC3 )
#define LAPACK_EXPORT_dtrevc F77_FUNC( dtrevc , DTREVC )
#define LAPACK_EXPORT_dtrexc F77_FUNC( dtrexc , DTREXC )
#define LAPACK_EXPORT_dtrrfs F77_FUNC( dtrrfs , DTRRFS )
#define LAPACK_EXPORT_dtrsen F77_FUNC( dtrsen , DTRSEN )
#define LAPACK_EXPORT_dtrsna F77_FUNC( dtrsna , DTRSNA )
#define LAPACK_EXPORT_dtrsyl F77_FUNC( dtrsyl , DTRSYL )
#define LAPACK_EXPORT_dtrti2 F77_FUNC( dtrti2 , DTRTI2 )
#define LAPACK_EXPORT_dtrtri F77_FUNC( dtrtri , DTRTRI )
#define LAPACK_EXPORT_dtrtrs F77_FUNC( dtrtrs , DTRTRS )
#define LAPACK_EXPORT_dtrttf F77_FUNC( dtrttf , DTRTTF )
#define LAPACK_EXPORT_dtrttp F77_FUNC( dtrttp , DTRTTP )
#define LAPACK_EXPORT_dtzrqf F77_FUNC( dtzrqf , DTZRQF )
#define LAPACK_EXPORT_dtzrzf F77_FUNC( dtzrzf , DTZRZF )
#define LAPACK_EXPORT_dzsum1 F77_FUNC( dzsum1 , DZSUM1 )
#define LAPACK_EXPORT_icmax1 F77_FUNC( icmax1 , ICMAX1 )
#define LAPACK_EXPORT_ilaclc F77_FUNC( ilaclc , ILACLC )
#define LAPACK_EXPORT_ilaclr F77_FUNC( ilaclr , ILACLR )
#define LAPACK_EXPORT_iladiag F77_FUNC( iladiag , ILADIAG )
#define LAPACK_EXPORT_iladlc F77_FUNC( iladlc , ILADLC )
#define LAPACK_EXPORT_iladlr F77_FUNC( iladlr , ILADLR )
#define LAPACK_EXPORT_ilaenv2stage F77_FUNC( ilaenv2stage , ILAENV2STAGE )
#define LAPACK_EXPORT_ilaprec F77_FUNC( ilaprec , ILAPREC )
#define LAPACK_EXPORT_ilaslc F77_FUNC( ilaslc , ILASLC )
#define LAPACK_EXPORT_ilaslr F77_FUNC( ilaslr , ILASLR )
#define LAPACK_EXPORT_ilatrans F77_FUNC( ilatrans , ILATRANS )
#define LAPACK_EXPORT_ilauplo F77_FUNC( ilauplo , ILAUPLO )
#define LAPACK_EXPORT_ilaver F77_FUNC( ilaver , ILAVER )
#define LAPACK_EXPORT_ilazlc F77_FUNC( ilazlc , ILAZLC )
#define LAPACK_EXPORT_ilazlr F77_FUNC( ilazlr , ILAZLR )
#define LAPACK_EXPORT_iparam2stage F77_FUNC( iparam2stage , IPARAM2STAGE )
#define LAPACK_EXPORT_izmax1 F77_FUNC( izmax1 , IZMAX1 )
#define LAPACK_EXPORT_sbbcsd F77_FUNC( sbbcsd , SBBCSD )
#define LAPACK_EXPORT_sbdsdc F77_FUNC( sbdsdc , SBDSDC )
#define LAPACK_EXPORT_sbdsqr F77_FUNC( sbdsqr , SBDSQR )
#define LAPACK_EXPORT_sbdsvdx F77_FUNC( sbdsvdx , SBDSVDX )
#define LAPACK_EXPORT_scombssq F77_FUNC( scombssq , SCOMBSSQ )
#define LAPACK_EXPORT_scsum1 F77_FUNC( scsum1 , SCSUM1 )
#define LAPACK_EXPORT_sdisna F77_FUNC( sdisna , SDISNA )
#define LAPACK_EXPORT_sgbbrd F77_FUNC( sgbbrd , SGBBRD )
#define LAPACK_EXPORT_sgbcon F77_FUNC( sgbcon , SGBCON )
#define LAPACK_EXPORT_sgbequb F77_FUNC( sgbequb , SGBEQUB )
#define LAPACK_EXPORT_sgbequ F77_FUNC( sgbequ , SGBEQU )
#define LAPACK_EXPORT_sgbrfs F77_FUNC( sgbrfs , SGBRFS )
#define LAPACK_EXPORT_sgbrfsx F77_FUNC( sgbrfsx , SGBRFSX )
#define LAPACK_EXPORT_sgbsv F77_FUNC( sgbsv , SGBSV )
#define LAPACK_EXPORT_sgbsvx F77_FUNC( sgbsvx , SGBSVX )
#define LAPACK_EXPORT_sgbsvxx F77_FUNC( sgbsvxx , SGBSVXX )
#define LAPACK_EXPORT_sgbtf2 F77_FUNC( sgbtf2 , SGBTF2 )
#define LAPACK_EXPORT_sgbtrf F77_FUNC( sgbtrf , SGBTRF )
#define LAPACK_EXPORT_sgbtrs F77_FUNC( sgbtrs , SGBTRS )
#define LAPACK_EXPORT_sgebak F77_FUNC( sgebak , SGEBAK )
#define LAPACK_EXPORT_sgebal F77_FUNC( sgebal , SGEBAL )
#define LAPACK_EXPORT_sgebd2 F77_FUNC( sgebd2 , SGEBD2 )
#define LAPACK_EXPORT_sgebrd F77_FUNC( sgebrd , SGEBRD )
#define LAPACK_EXPORT_sgecon F77_FUNC( sgecon , SGECON )
#define LAPACK_EXPORT_sgeequb F77_FUNC( sgeequb , SGEEQUB )
#define LAPACK_EXPORT_sgeequ F77_FUNC( sgeequ , SGEEQU )
#define LAPACK_EXPORT_sgees F77_FUNC( sgees , SGEES )
#define LAPACK_EXPORT_sgeesx F77_FUNC( sgeesx , SGEESX )
#define LAPACK_EXPORT_sgeev F77_FUNC( sgeev , SGEEV )
#define LAPACK_EXPORT_sgeevx F77_FUNC( sgeevx , SGEEVX )
#define LAPACK_EXPORT_sgegs F77_FUNC( sgegs , SGEGS )
#define LAPACK_EXPORT_sgegv F77_FUNC( sgegv , SGEGV )
#define LAPACK_EXPORT_sgehd2 F77_FUNC( sgehd2 , SGEHD2 )
#define LAPACK_EXPORT_sgehrd F77_FUNC( sgehrd , SGEHRD )
#define LAPACK_EXPORT_sgejsv F77_FUNC( sgejsv , SGEJSV )
#define LAPACK_EXPORT_sgelq2 F77_FUNC( sgelq2 , SGELQ2 )
#define LAPACK_EXPORT_sgelq F77_FUNC( sgelq , SGELQ )
#define LAPACK_EXPORT_sgelqf F77_FUNC( sgelqf , SGELQF )
#define LAPACK_EXPORT_sgelqt3 F77_FUNC( sgelqt3 , SGELQT3 )
#define LAPACK_EXPORT_sgelqt F77_FUNC( sgelqt , SGELQT )
#define LAPACK_EXPORT_sgelsd F77_FUNC( sgelsd , SGELSD )
#define LAPACK_EXPORT_sgels F77_FUNC( sgels , SGELS )
#define LAPACK_EXPORT_sgelss F77_FUNC( sgelss , SGELSS )
#define LAPACK_EXPORT_sgelsx F77_FUNC( sgelsx , SGELSX )
#define LAPACK_EXPORT_sgelsy F77_FUNC( sgelsy , SGELSY )
#define LAPACK_EXPORT_sgemlq F77_FUNC( sgemlq , SGEMLQ )
#define LAPACK_EXPORT_sgemlqt F77_FUNC( sgemlqt , SGEMLQT )
#define LAPACK_EXPORT_sgemqr F77_FUNC( sgemqr , SGEMQR )
#define LAPACK_EXPORT_sgemqrt F77_FUNC( sgemqrt , SGEMQRT )
#define LAPACK_EXPORT_sgeql2 F77_FUNC( sgeql2 , SGEQL2 )
#define LAPACK_EXPORT_sgeqlf F77_FUNC( sgeqlf , SGEQLF )
#define LAPACK_EXPORT_sgeqp3 F77_FUNC( sgeqp3 , SGEQP3 )
#define LAPACK_EXPORT_sgeqpf F77_FUNC( sgeqpf , SGEQPF )
#define LAPACK_EXPORT_sgeqr2 F77_FUNC( sgeqr2 , SGEQR2 )
#define LAPACK_EXPORT_sgeqr2p F77_FUNC( sgeqr2p , SGEQR2P )
#define LAPACK_EXPORT_sgeqr F77_FUNC( sgeqr , SGEQR )
#define LAPACK_EXPORT_sgeqrf F77_FUNC( sgeqrf , SGEQRF )
#define LAPACK_EXPORT_sgeqrfp F77_FUNC( sgeqrfp , SGEQRFP )
#define LAPACK_EXPORT_sgeqrt2 F77_FUNC( sgeqrt2 , SGEQRT2 )
#define LAPACK_EXPORT_sgeqrt3 F77_FUNC( sgeqrt3 , SGEQRT3 )
#define LAPACK_EXPORT_sgeqrt F77_FUNC( sgeqrt , SGEQRT )
#define LAPACK_EXPORT_sgerfs F77_FUNC( sgerfs , SGERFS )
#define LAPACK_EXPORT_sgerfsx F77_FUNC( sgerfsx , SGERFSX )
#define LAPACK_EXPORT_sgerq2 F77_FUNC( sgerq2 , SGERQ2 )
#define LAPACK_EXPORT_sgerqf F77_FUNC( sgerqf , SGERQF )
#define LAPACK_EXPORT_sgesc2 F77_FUNC( sgesc2 , SGESC2 )
#define LAPACK_EXPORT_sgesdd F77_FUNC( sgesdd , SGESDD )
#define LAPACK_EXPORT_sgesvd F77_FUNC( sgesvd , SGESVD )
#define LAPACK_EXPORT_sgesvdq F77_FUNC( sgesvdq , SGESVDQ )
#define LAPACK_EXPORT_sgesvdx F77_FUNC( sgesvdx , SGESVDX )
#define LAPACK_EXPORT_sgesv F77_FUNC( sgesv , SGESV )
#define LAPACK_EXPORT_sgesvj F77_FUNC( sgesvj , SGESVJ )
#define LAPACK_EXPORT_sgesvx F77_FUNC( sgesvx , SGESVX )
#define LAPACK_EXPORT_sgesvxx F77_FUNC( sgesvxx , SGESVXX )
#define LAPACK_EXPORT_sgetc2 F77_FUNC( sgetc2 , SGETC2 )
#define LAPACK_EXPORT_sgetf2 F77_FUNC( sgetf2 , SGETF2 )
#define LAPACK_EXPORT_sgetrf2 F77_FUNC( sgetrf2 , SGETRF2 )
#define LAPACK_EXPORT_sgetrf F77_FUNC( sgetrf , SGETRF )
#define LAPACK_EXPORT_sgetri F77_FUNC( sgetri , SGETRI )
#define LAPACK_EXPORT_sgetrs F77_FUNC( sgetrs , SGETRS )
#define LAPACK_EXPORT_sgetsls F77_FUNC( sgetsls , SGETSLS )
#define LAPACK_EXPORT_sggbak F77_FUNC( sggbak , SGGBAK )
#define LAPACK_EXPORT_sggbal F77_FUNC( sggbal , SGGBAL )
#define LAPACK_EXPORT_sgges3 F77_FUNC( sgges3 , SGGES3 )
#define LAPACK_EXPORT_sgges F77_FUNC( sgges , SGGES )
#define LAPACK_EXPORT_sggesx F77_FUNC( sggesx , SGGESX )
#define LAPACK_EXPORT_sggev3 F77_FUNC( sggev3 , SGGEV3 )
#define LAPACK_EXPORT_sggev F77_FUNC( sggev , SGGEV )
#define LAPACK_EXPORT_sggevx F77_FUNC( sggevx , SGGEVX )
#define LAPACK_EXPORT_sggglm F77_FUNC( sggglm , SGGGLM )
#define LAPACK_EXPORT_sgghd3 F77_FUNC( sgghd3 , SGGHD3 )
#define LAPACK_EXPORT_sgghrd F77_FUNC( sgghrd , SGGHRD )
#define LAPACK_EXPORT_sgglse F77_FUNC( sgglse , SGGLSE )
#define LAPACK_EXPORT_sggqrf F77_FUNC( sggqrf , SGGQRF )
#define LAPACK_EXPORT_sggrqf F77_FUNC( sggrqf , SGGRQF )
#define LAPACK_EXPORT_sggsvd3 F77_FUNC( sggsvd3 , SGGSVD3 )
#define LAPACK_EXPORT_sggsvd F77_FUNC( sggsvd , SGGSVD )
#define LAPACK_EXPORT_sggsvp3 F77_FUNC( sggsvp3 , SGGSVP3 )
#define LAPACK_EXPORT_sggsvp F77_FUNC( sggsvp , SGGSVP )
#define LAPACK_EXPORT_sgsvj0 F77_FUNC( sgsvj0 , SGSVJ0 )
#define LAPACK_EXPORT_sgsvj1 F77_FUNC( sgsvj1 , SGSVJ1 )
#define LAPACK_EXPORT_sgtcon F77_FUNC( sgtcon , SGTCON )
#define LAPACK_EXPORT_sgtrfs F77_FUNC( sgtrfs , SGTRFS )
#define LAPACK_EXPORT_sgtsv F77_FUNC( sgtsv , SGTSV )
#define LAPACK_EXPORT_sgtsvx F77_FUNC( sgtsvx , SGTSVX )
#define LAPACK_EXPORT_sgttrf F77_FUNC( sgttrf , SGTTRF )
#define LAPACK_EXPORT_sgttrs F77_FUNC( sgttrs , SGTTRS )
#define LAPACK_EXPORT_sgtts2 F77_FUNC( sgtts2 , SGTTS2 )
#define LAPACK_EXPORT_shgeqz F77_FUNC( shgeqz , SHGEQZ )
#define LAPACK_EXPORT_shsein F77_FUNC( shsein , SHSEIN )
#define LAPACK_EXPORT_shseqr F77_FUNC( shseqr , SHSEQR )
#define LAPACK_EXPORT_sisnan F77_FUNC( sisnan , SISNAN )
#define LAPACK_EXPORT_slabad F77_FUNC( slabad , SLABAD )
#define LAPACK_EXPORT_slabrd F77_FUNC( slabrd , SLABRD )
#define LAPACK_EXPORT_slacn2 F77_FUNC( slacn2 , SLACN2 )
#define LAPACK_EXPORT_slacon F77_FUNC( slacon , SLACON )
#define LAPACK_EXPORT_slacpy F77_FUNC( slacpy , SLACPY )
#define LAPACK_EXPORT_sladiv F77_FUNC( sladiv , SLADIV )
#define LAPACK_EXPORT_slae2 F77_FUNC( slae2 , SLAE2 )
#define LAPACK_EXPORT_slaebz F77_FUNC( slaebz , SLAEBZ )
#define LAPACK_EXPORT_slaed0 F77_FUNC( slaed0 , SLAED0 )
#define LAPACK_EXPORT_slaed1 F77_FUNC( slaed1 , SLAED1 )
#define LAPACK_EXPORT_slaed2 F77_FUNC( slaed2 , SLAED2 )
#define LAPACK_EXPORT_slaed3 F77_FUNC( slaed3 , SLAED3 )
#define LAPACK_EXPORT_slaed4 F77_FUNC( slaed4 , SLAED4 )
#define LAPACK_EXPORT_slaed5 F77_FUNC( slaed5 , SLAED5 )
#define LAPACK_EXPORT_slaed6 F77_FUNC( slaed6 , SLAED6 )
#define LAPACK_EXPORT_slaed7 F77_FUNC( slaed7 , SLAED7 )
#define LAPACK_EXPORT_slaed8 F77_FUNC( slaed8 , SLAED8 )
#define LAPACK_EXPORT_slaed9 F77_FUNC( slaed9 , SLAED9 )
#define LAPACK_EXPORT_slaeda F77_FUNC( slaeda , SLAEDA )
#define LAPACK_EXPORT_slaein F77_FUNC( slaein , SLAEIN )
#define LAPACK_EXPORT_slaev2 F77_FUNC( slaev2 , SLAEV2 )
#define LAPACK_EXPORT_slaexc F77_FUNC( slaexc , SLAEXC )
#define LAPACK_EXPORT_slag2d F77_FUNC( slag2d , SLAG2D )
#define LAPACK_EXPORT_slag2 F77_FUNC( slag2 , SLAG2 )
#define LAPACK_EXPORT_sla_gbamv F77_FUNC( sla_gbamv , SLA_GBAMV )
#define LAPACK_EXPORT_sla_gbrcond F77_FUNC( sla_gbrcond , SLA_GBRCOND )
#define LAPACK_EXPORT_sla_gbrfsx_extended F77_FUNC( sla_gbrfsx_extended , SLA_GBRFSX_EXTENDED )
#define LAPACK_EXPORT_sla_gbrpvgrw F77_FUNC( sla_gbrpvgrw , SLA_GBRPVGRW )
#define LAPACK_EXPORT_sla_geamv F77_FUNC( sla_geamv , SLA_GEAMV )
#define LAPACK_EXPORT_sla_gercond F77_FUNC( sla_gercond , SLA_GERCOND )
#define LAPACK_EXPORT_sla_gerfsx_extended F77_FUNC( sla_gerfsx_extended , SLA_GERFSX_EXTENDED )
#define LAPACK_EXPORT_sla_gerpvgrw F77_FUNC( sla_gerpvgrw , SLA_GERPVGRW )
#define LAPACK_EXPORT_slags2 F77_FUNC( slags2 , SLAGS2 )
#define LAPACK_EXPORT_slagtf F77_FUNC( slagtf , SLAGTF )
#define LAPACK_EXPORT_slagtm F77_FUNC( slagtm , SLAGTM )
#define LAPACK_EXPORT_slagts F77_FUNC( slagts , SLAGTS )
#define LAPACK_EXPORT_slagv2 F77_FUNC( slagv2 , SLAGV2 )
#define LAPACK_EXPORT_slahqr F77_FUNC( slahqr , SLAHQR )
#define LAPACK_EXPORT_slahr2 F77_FUNC( slahr2 , SLAHR2 )
#define LAPACK_EXPORT_slahrd F77_FUNC( slahrd , SLAHRD )
#define LAPACK_EXPORT_slaic1 F77_FUNC( slaic1 , SLAIC1 )
#define LAPACK_EXPORT_slaisnan F77_FUNC( slaisnan , SLAISNAN )
#define LAPACK_EXPORT_sla_lin_berr F77_FUNC( sla_lin_berr , SLA_LIN_BERR )
#define LAPACK_EXPORT_slaln2 F77_FUNC( slaln2 , SLALN2 )
#define LAPACK_EXPORT_slals0 F77_FUNC( slals0 , SLALS0 )
#define LAPACK_EXPORT_slalsa F77_FUNC( slalsa , SLALSA )
#define LAPACK_EXPORT_slalsd F77_FUNC( slalsd , SLALSD )
#define LAPACK_EXPORT_slamrg F77_FUNC( slamrg , SLAMRG )
#define LAPACK_EXPORT_slamswlq F77_FUNC( slamswlq , SLAMSWLQ )
#define LAPACK_EXPORT_slamtsqr F77_FUNC( slamtsqr , SLAMTSQR )
#define LAPACK_EXPORT_slaneg F77_FUNC( slaneg , SLANEG )
#define LAPACK_EXPORT_slangb F77_FUNC( slangb , SLANGB )
#define LAPACK_EXPORT_slange F77_FUNC( slange , SLANGE )
#define LAPACK_EXPORT_slangt F77_FUNC( slangt , SLANGT )
#define LAPACK_EXPORT_slanhs F77_FUNC( slanhs , SLANHS )
#define LAPACK_EXPORT_slansb F77_FUNC( slansb , SLANSB )
#define LAPACK_EXPORT_slansf F77_FUNC( slansf , SLANSF )
#define LAPACK_EXPORT_slansp F77_FUNC( slansp , SLANSP )
#define LAPACK_EXPORT_slanst F77_FUNC( slanst , SLANST )
#define LAPACK_EXPORT_slansy F77_FUNC( slansy , SLANSY )
#define LAPACK_EXPORT_slantb F77_FUNC( slantb , SLANTB )
#define LAPACK_EXPORT_slantp F77_FUNC( slantp , SLANTP )
#define LAPACK_EXPORT_slantr F77_FUNC( slantr , SLANTR )
#define LAPACK_EXPORT_slanv2 F77_FUNC( slanv2 , SLANV2 )
#define LAPACK_EXPORT_slaorhr_col_getrfnp2 F77_FUNC( slaorhr_col_getrfnp2 , SLAORHR_COL_GETRFNP2 )
#define LAPACK_EXPORT_slaorhr_col_getrfnp F77_FUNC( slaorhr_col_getrfnp , SLAORHR_COL_GETRFNP )
#define LAPACK_EXPORT_slapll F77_FUNC( slapll , SLAPLL )
#define LAPACK_EXPORT_slapmr F77_FUNC( slapmr , SLAPMR )
#define LAPACK_EXPORT_slapmt F77_FUNC( slapmt , SLAPMT )
#define LAPACK_EXPORT_sla_porcond F77_FUNC( sla_porcond , SLA_PORCOND )
#define LAPACK_EXPORT_sla_porfsx_extended F77_FUNC( sla_porfsx_extended , SLA_PORFSX_EXTENDED )
#define LAPACK_EXPORT_sla_porpvgrw F77_FUNC( sla_porpvgrw , SLA_PORPVGRW )
#define LAPACK_EXPORT_slapy2 F77_FUNC( slapy2 , SLAPY2 )
#define LAPACK_EXPORT_slapy3 F77_FUNC( slapy3 , SLAPY3 )
#define LAPACK_EXPORT_slaqgb F77_FUNC( slaqgb , SLAQGB )
#define LAPACK_EXPORT_slaqge F77_FUNC( slaqge , SLAQGE )
#define LAPACK_EXPORT_slaqp2 F77_FUNC( slaqp2 , SLAQP2 )
#define LAPACK_EXPORT_slaqps F77_FUNC( slaqps , SLAQPS )
#define LAPACK_EXPORT_slaqr0 F77_FUNC( slaqr0 , SLAQR0 )
#define LAPACK_EXPORT_slaqr1 F77_FUNC( slaqr1 , SLAQR1 )
#define LAPACK_EXPORT_slaqr2 F77_FUNC( slaqr2 , SLAQR2 )
#define LAPACK_EXPORT_slaqr3 F77_FUNC( slaqr3 , SLAQR3 )
#define LAPACK_EXPORT_slaqr4 F77_FUNC( slaqr4 , SLAQR4 )
#define LAPACK_EXPORT_slaqr5 F77_FUNC( slaqr5 , SLAQR5 )
#define LAPACK_EXPORT_slaqsb F77_FUNC( slaqsb , SLAQSB )
#define LAPACK_EXPORT_slaqsp F77_FUNC( slaqsp , SLAQSP )
#define LAPACK_EXPORT_slaqsy F77_FUNC( slaqsy , SLAQSY )
#define LAPACK_EXPORT_slaqtr F77_FUNC( slaqtr , SLAQTR )
#define LAPACK_EXPORT_slar1v F77_FUNC( slar1v , SLAR1V )
#define LAPACK_EXPORT_slar2v F77_FUNC( slar2v , SLAR2V )
#define LAPACK_EXPORT_slarfb F77_FUNC( slarfb , SLARFB )
#define LAPACK_EXPORT_slarf F77_FUNC( slarf , SLARF )
#define LAPACK_EXPORT_slarfg F77_FUNC( slarfg , SLARFG )
#define LAPACK_EXPORT_slarfgp F77_FUNC( slarfgp , SLARFGP )
#define LAPACK_EXPORT_slarft F77_FUNC( slarft , SLARFT )
#define LAPACK_EXPORT_slarfx F77_FUNC( slarfx , SLARFX )
#define LAPACK_EXPORT_slarfy F77_FUNC( slarfy , SLARFY )
#define LAPACK_EXPORT_slargv F77_FUNC( slargv , SLARGV )
#define LAPACK_EXPORT_slarnv F77_FUNC( slarnv , SLARNV )
#define LAPACK_EXPORT_slarra F77_FUNC( slarra , SLARRA )
#define LAPACK_EXPORT_slarrb F77_FUNC( slarrb , SLARRB )
#define LAPACK_EXPORT_slarrc F77_FUNC( slarrc , SLARRC )
#define LAPACK_EXPORT_slarrd F77_FUNC( slarrd , SLARRD )
#define LAPACK_EXPORT_slarre F77_FUNC( slarre , SLARRE )
#define LAPACK_EXPORT_slarrf F77_FUNC( slarrf , SLARRF )
#define LAPACK_EXPORT_slarrj F77_FUNC( slarrj , SLARRJ )
#define LAPACK_EXPORT_slarrk F77_FUNC( slarrk , SLARRK )
#define LAPACK_EXPORT_slarrr F77_FUNC( slarrr , SLARRR )
#define LAPACK_EXPORT_slarrv F77_FUNC( slarrv , SLARRV )
#define LAPACK_EXPORT_slarscl2 F77_FUNC( slarscl2 , SLARSCL2 )
#define LAPACK_EXPORT_slartg F77_FUNC( slartg , SLARTG )
#define LAPACK_EXPORT_slartgp F77_FUNC( slartgp , SLARTGP )
#define LAPACK_EXPORT_slartgs F77_FUNC( slartgs , SLARTGS )
#define LAPACK_EXPORT_slartv F77_FUNC( slartv , SLARTV )
#define LAPACK_EXPORT_slaruv F77_FUNC( slaruv , SLARUV )
#define LAPACK_EXPORT_slarzb F77_FUNC( slarzb , SLARZB )
#define LAPACK_EXPORT_slarz F77_FUNC( slarz , SLARZ )
#define LAPACK_EXPORT_slarzt F77_FUNC( slarzt , SLARZT )
#define LAPACK_EXPORT_slas2 F77_FUNC( slas2 , SLAS2 )
#define LAPACK_EXPORT_slascl2 F77_FUNC( slascl2 , SLASCL2 )
#define LAPACK_EXPORT_slascl F77_FUNC( slascl , SLASCL )
#define LAPACK_EXPORT_slasd0 F77_FUNC( slasd0 , SLASD0 )
#define LAPACK_EXPORT_slasd1 F77_FUNC( slasd1 , SLASD1 )
#define LAPACK_EXPORT_slasd2 F77_FUNC( slasd2 , SLASD2 )
#define LAPACK_EXPORT_slasd3 F77_FUNC( slasd3 , SLASD3 )
#define LAPACK_EXPORT_slasd4 F77_FUNC( slasd4 , SLASD4 )
#define LAPACK_EXPORT_slasd5 F77_FUNC( slasd5 , SLASD5 )
#define LAPACK_EXPORT_slasd6 F77_FUNC( slasd6 , SLASD6 )
#define LAPACK_EXPORT_slasd7 F77_FUNC( slasd7 , SLASD7 )
#define LAPACK_EXPORT_slasd8 F77_FUNC( slasd8 , SLASD8 )
#define LAPACK_EXPORT_slasda F77_FUNC( slasda , SLASDA )
#define LAPACK_EXPORT_slasdq F77_FUNC( slasdq , SLASDQ )
#define LAPACK_EXPORT_slasdt F77_FUNC( slasdt , SLASDT )
#define LAPACK_EXPORT_slaset F77_FUNC( slaset , SLASET )
#define LAPACK_EXPORT_slasq1 F77_FUNC( slasq1 , SLASQ1 )
#define LAPACK_EXPORT_slasq2 F77_FUNC( slasq2 , SLASQ2 )
#define LAPACK_EXPORT_slasq3 F77_FUNC( slasq3 , SLASQ3 )
#define LAPACK_EXPORT_slasq4 F77_FUNC( slasq4 , SLASQ4 )
#define LAPACK_EXPORT_slasq5 F77_FUNC( slasq5 , SLASQ5 )
#define LAPACK_EXPORT_slasq6 F77_FUNC( slasq6 , SLASQ6 )
#define LAPACK_EXPORT_slasr F77_FUNC( slasr , SLASR )
#define LAPACK_EXPORT_slasrt F77_FUNC( slasrt , SLASRT )
#define LAPACK_EXPORT_slassq F77_FUNC( slassq , SLASSQ )
#define LAPACK_EXPORT_slasv2 F77_FUNC( slasv2 , SLASV2 )
#define LAPACK_EXPORT_slaswlq F77_FUNC( slaswlq , SLASWLQ )
#define LAPACK_EXPORT_slaswp F77_FUNC( slaswp , SLASWP )
#define LAPACK_EXPORT_slasy2 F77_FUNC( slasy2 , SLASY2 )
#define LAPACK_EXPORT_sla_syamv F77_FUNC( sla_syamv , SLA_SYAMV )
#define LAPACK_EXPORT_slasyf_aa F77_FUNC( slasyf_aa , SLASYF_AA )
#define LAPACK_EXPORT_slasyf F77_FUNC( slasyf , SLASYF )
#define LAPACK_EXPORT_slasyf_rk F77_FUNC( slasyf_rk , SLASYF_RK )
#define LAPACK_EXPORT_slasyf_rook F77_FUNC( slasyf_rook , SLASYF_ROOK )
#define LAPACK_EXPORT_sla_syrcond F77_FUNC( sla_syrcond , SLA_SYRCOND )
#define LAPACK_EXPORT_sla_syrfsx_extended F77_FUNC( sla_syrfsx_extended , SLA_SYRFSX_EXTENDED )
#define LAPACK_EXPORT_sla_syrpvgrw F77_FUNC( sla_syrpvgrw , SLA_SYRPVGRW )
#define LAPACK_EXPORT_slatbs F77_FUNC( slatbs , SLATBS )
#define LAPACK_EXPORT_slatdf F77_FUNC( slatdf , SLATDF )
#define LAPACK_EXPORT_slatps F77_FUNC( slatps , SLATPS )
#define LAPACK_EXPORT_slatrd F77_FUNC( slatrd , SLATRD )
#define LAPACK_EXPORT_slatrs F77_FUNC( slatrs , SLATRS )
#define LAPACK_EXPORT_slatrz F77_FUNC( slatrz , SLATRZ )
#define LAPACK_EXPORT_slatsqr F77_FUNC( slatsqr , SLATSQR )
#define LAPACK_EXPORT_slatzm F77_FUNC( slatzm , SLATZM )
#define LAPACK_EXPORT_slauu2 F77_FUNC( slauu2 , SLAUU2 )
#define LAPACK_EXPORT_slauum F77_FUNC( slauum , SLAUUM )
#define LAPACK_EXPORT_sla_wwaddw F77_FUNC( sla_wwaddw , SLA_WWADDW )
#define LAPACK_EXPORT_sopmtr F77_FUNC( sopmtr , SOPMTR )
#define LAPACK_EXPORT_sorbdb1 F77_FUNC( sorbdb1 , SORBDB1 )
#define LAPACK_EXPORT_sorbdb2 F77_FUNC( sorbdb2 , SORBDB2 )
#define LAPACK_EXPORT_sorbdb3 F77_FUNC( sorbdb3 , SORBDB3 )
#define LAPACK_EXPORT_sorbdb4 F77_FUNC( sorbdb4 , SORBDB4 )
#define LAPACK_EXPORT_sorbdb5 F77_FUNC( sorbdb5 , SORBDB5 )
#define LAPACK_EXPORT_sorbdb6 F77_FUNC( sorbdb6 , SORBDB6 )
#define LAPACK_EXPORT_sorbdb F77_FUNC( sorbdb , SORBDB )
#define LAPACK_EXPORT_sorg2l F77_FUNC( sorg2l , SORG2L )
#define LAPACK_EXPORT_sorg2r F77_FUNC( sorg2r , SORG2R )
#define LAPACK_EXPORT_sorgbr F77_FUNC( sorgbr , SORGBR )
#define LAPACK_EXPORT_sorgl2 F77_FUNC( sorgl2 , SORGL2 )
#define LAPACK_EXPORT_sorglq F77_FUNC( sorglq , SORGLQ )
#define LAPACK_EXPORT_sorgql F77_FUNC( sorgql , SORGQL )
#define LAPACK_EXPORT_sorgqr F77_FUNC( sorgqr , SORGQR )
#define LAPACK_EXPORT_sorgr2 F77_FUNC( sorgr2 , SORGR2 )
#define LAPACK_EXPORT_sorgrq F77_FUNC( sorgrq , SORGRQ )
#define LAPACK_EXPORT_sorgtr F77_FUNC( sorgtr , SORGTR )
#define LAPACK_EXPORT_sorgtsqr F77_FUNC( sorgtsqr , SORGTSQR )
#define LAPACK_EXPORT_sorhr_col F77_FUNC( sorhr_col , SORHR_COL )
#define LAPACK_EXPORT_sorm22 F77_FUNC( sorm22 , SORM22 )
#define LAPACK_EXPORT_sorm2l F77_FUNC( sorm2l , SORM2L )
#define LAPACK_EXPORT_sorm2r F77_FUNC( sorm2r , SORM2R )
#define LAPACK_EXPORT_sormbr F77_FUNC( sormbr , SORMBR )
#define LAPACK_EXPORT_sorml2 F77_FUNC( sorml2 , SORML2 )
#define LAPACK_EXPORT_sormlq F77_FUNC( sormlq , SORMLQ )
#define LAPACK_EXPORT_sormql F77_FUNC( sormql , SORMQL )
#define LAPACK_EXPORT_sormqr F77_FUNC( sormqr , SORMQR )
#define LAPACK_EXPORT_sormr2 F77_FUNC( sormr2 , SORMR2 )
#define LAPACK_EXPORT_sormr3 F77_FUNC( sormr3 , SORMR3 )
#define LAPACK_EXPORT_sormrq F77_FUNC( sormrq , SORMRQ )
#define LAPACK_EXPORT_sormrz F77_FUNC( sormrz , SORMRZ )
#define LAPACK_EXPORT_sormtr F77_FUNC( sormtr , SORMTR )
#define LAPACK_EXPORT_spbcon F77_FUNC( spbcon , SPBCON )
#define LAPACK_EXPORT_spbequ F77_FUNC( spbequ , SPBEQU )
#define LAPACK_EXPORT_spbrfs F77_FUNC( spbrfs , SPBRFS )
#define LAPACK_EXPORT_spbstf F77_FUNC( spbstf , SPBSTF )
#define LAPACK_EXPORT_spbsv F77_FUNC( spbsv , SPBSV )
#define LAPACK_EXPORT_spbsvx F77_FUNC( spbsvx , SPBSVX )
#define LAPACK_EXPORT_spbtf2 F77_FUNC( spbtf2 , SPBTF2 )
#define LAPACK_EXPORT_spbtrf F77_FUNC( spbtrf , SPBTRF )
#define LAPACK_EXPORT_spbtrs F77_FUNC( spbtrs , SPBTRS )
#define LAPACK_EXPORT_spftrf F77_FUNC( spftrf , SPFTRF )
#define LAPACK_EXPORT_spftri F77_FUNC( spftri , SPFTRI )
#define LAPACK_EXPORT_spftrs F77_FUNC( spftrs , SPFTRS )
#define LAPACK_EXPORT_spocon F77_FUNC( spocon , SPOCON )
#define LAPACK_EXPORT_spoequb F77_FUNC( spoequb , SPOEQUB )
#define LAPACK_EXPORT_spoequ F77_FUNC( spoequ , SPOEQU )
#define LAPACK_EXPORT_sporfs F77_FUNC( sporfs , SPORFS )
#define LAPACK_EXPORT_sporfsx F77_FUNC( sporfsx , SPORFSX )
#define LAPACK_EXPORT_sposv F77_FUNC( sposv , SPOSV )
#define LAPACK_EXPORT_sposvx F77_FUNC( sposvx , SPOSVX )
#define LAPACK_EXPORT_sposvxx F77_FUNC( sposvxx , SPOSVXX )
#define LAPACK_EXPORT_spotf2 F77_FUNC( spotf2 , SPOTF2 )
#define LAPACK_EXPORT_spotrf2 F77_FUNC( spotrf2 , SPOTRF2 )
#define LAPACK_EXPORT_spotrf F77_FUNC( spotrf , SPOTRF )
#define LAPACK_EXPORT_spotri F77_FUNC( spotri , SPOTRI )
#define LAPACK_EXPORT_spotrs F77_FUNC( spotrs , SPOTRS )
#define LAPACK_EXPORT_sppcon F77_FUNC( sppcon , SPPCON )
#define LAPACK_EXPORT_sppequ F77_FUNC( sppequ , SPPEQU )
#define LAPACK_EXPORT_spprfs F77_FUNC( spprfs , SPPRFS )
#define LAPACK_EXPORT_sppsv F77_FUNC( sppsv , SPPSV )
#define LAPACK_EXPORT_sppsvx F77_FUNC( sppsvx , SPPSVX )
#define LAPACK_EXPORT_spptrf F77_FUNC( spptrf , SPPTRF )
#define LAPACK_EXPORT_spptri F77_FUNC( spptri , SPPTRI )
#define LAPACK_EXPORT_spptrs F77_FUNC( spptrs , SPPTRS )
#define LAPACK_EXPORT_spstf2 F77_FUNC( spstf2 , SPSTF2 )
#define LAPACK_EXPORT_spstrf F77_FUNC( spstrf , SPSTRF )
#define LAPACK_EXPORT_sptcon F77_FUNC( sptcon , SPTCON )
#define LAPACK_EXPORT_spteqr F77_FUNC( spteqr , SPTEQR )
#define LAPACK_EXPORT_sptrfs F77_FUNC( sptrfs , SPTRFS )
#define LAPACK_EXPORT_sptsv F77_FUNC( sptsv , SPTSV )
#define LAPACK_EXPORT_sptsvx F77_FUNC( sptsvx , SPTSVX )
#define LAPACK_EXPORT_spttrf F77_FUNC( spttrf , SPTTRF )
#define LAPACK_EXPORT_spttrs F77_FUNC( spttrs , SPTTRS )
#define LAPACK_EXPORT_sptts2 F77_FUNC( sptts2 , SPTTS2 )
#define LAPACK_EXPORT_srscl F77_FUNC( srscl , SRSCL )
#define LAPACK_EXPORT_ssb2st_kernels F77_FUNC( ssb2st_kernels , SSB2ST_KERNELS )
#define LAPACK_EXPORT_ssbev_2stage F77_FUNC( ssbev_2stage , SSBEV_2STAGE )
#define LAPACK_EXPORT_ssbevd_2stage F77_FUNC( ssbevd_2stage , SSBEVD_2STAGE )
#define LAPACK_EXPORT_ssbevd F77_FUNC( ssbevd , SSBEVD )
#define LAPACK_EXPORT_ssbev F77_FUNC( ssbev , SSBEV )
#define LAPACK_EXPORT_ssbevx_2stage F77_FUNC( ssbevx_2stage , SSBEVX_2STAGE )
#define LAPACK_EXPORT_ssbevx F77_FUNC( ssbevx , SSBEVX )
#define LAPACK_EXPORT_ssbgst F77_FUNC( ssbgst , SSBGST )
#define LAPACK_EXPORT_ssbgvd F77_FUNC( ssbgvd , SSBGVD )
#define LAPACK_EXPORT_ssbgv F77_FUNC( ssbgv , SSBGV )
#define LAPACK_EXPORT_ssbgvx F77_FUNC( ssbgvx , SSBGVX )
#define LAPACK_EXPORT_ssbtrd F77_FUNC( ssbtrd , SSBTRD )
#define LAPACK_EXPORT_ssfrk F77_FUNC( ssfrk , SSFRK )
#define LAPACK_EXPORT_sspcon F77_FUNC( sspcon , SSPCON )
#define LAPACK_EXPORT_sspevd F77_FUNC( sspevd , SSPEVD )
#define LAPACK_EXPORT_sspev F77_FUNC( sspev , SSPEV )
#define LAPACK_EXPORT_sspevx F77_FUNC( sspevx , SSPEVX )
#define LAPACK_EXPORT_sspgst F77_FUNC( sspgst , SSPGST )
#define LAPACK_EXPORT_sspgvd F77_FUNC( sspgvd , SSPGVD )
#define LAPACK_EXPORT_sspgv F77_FUNC( sspgv , SSPGV )
#define LAPACK_EXPORT_sspgvx F77_FUNC( sspgvx , SSPGVX )
#define LAPACK_EXPORT_ssprfs F77_FUNC( ssprfs , SSPRFS )
#define LAPACK_EXPORT_sspsv F77_FUNC( sspsv , SSPSV )
#define LAPACK_EXPORT_sspsvx F77_FUNC( sspsvx , SSPSVX )
#define LAPACK_EXPORT_ssptrd F77_FUNC( ssptrd , SSPTRD )
#define LAPACK_EXPORT_ssptrf F77_FUNC( ssptrf , SSPTRF )
#define LAPACK_EXPORT_ssptri F77_FUNC( ssptri , SSPTRI )
#define LAPACK_EXPORT_ssptrs F77_FUNC( ssptrs , SSPTRS )
#define LAPACK_EXPORT_sstebz F77_FUNC( sstebz , SSTEBZ )
#define LAPACK_EXPORT_sstedc F77_FUNC( sstedc , SSTEDC )
#define LAPACK_EXPORT_sstegr F77_FUNC( sstegr , SSTEGR )
#define LAPACK_EXPORT_sstein F77_FUNC( sstein , SSTEIN )
#define LAPACK_EXPORT_sstemr F77_FUNC( sstemr , SSTEMR )
#define LAPACK_EXPORT_ssteqr F77_FUNC( ssteqr , SSTEQR )
#define LAPACK_EXPORT_ssterf F77_FUNC( ssterf , SSTERF )
#define LAPACK_EXPORT_sstevd F77_FUNC( sstevd , SSTEVD )
#define LAPACK_EXPORT_sstev F77_FUNC( sstev , SSTEV )
#define LAPACK_EXPORT_sstevr F77_FUNC( sstevr , SSTEVR )
#define LAPACK_EXPORT_sstevx F77_FUNC( sstevx , SSTEVX )
#define LAPACK_EXPORT_ssycon_3 F77_FUNC( ssycon_3 , SSYCON_3 )
#define LAPACK_EXPORT_ssycon F77_FUNC( ssycon , SSYCON )
#define LAPACK_EXPORT_ssycon_rook F77_FUNC( ssycon_rook , SSYCON_ROOK )
#define LAPACK_EXPORT_ssyconv F77_FUNC( ssyconv , SSYCONV )
#define LAPACK_EXPORT_ssyconvf F77_FUNC( ssyconvf , SSYCONVF )
#define LAPACK_EXPORT_ssyconvf_rook F77_FUNC( ssyconvf_rook , SSYCONVF_ROOK )
#define LAPACK_EXPORT_ssyequb F77_FUNC( ssyequb , SSYEQUB )
#define LAPACK_EXPORT_ssyev_2stage F77_FUNC( ssyev_2stage , SSYEV_2STAGE )
#define LAPACK_EXPORT_ssyevd_2stage F77_FUNC( ssyevd_2stage , SSYEVD_2STAGE )
#define LAPACK_EXPORT_ssyevd F77_FUNC( ssyevd , SSYEVD )
#define LAPACK_EXPORT_ssyev  F77_FUNC( ssyev  , SSYEV  )
#define LAPACK_EXPORT_ssyev F77_FUNC( ssyev , SSYEV )
#define LAPACK_EXPORT_ssyevr_2stage F77_FUNC( ssyevr_2stage , SSYEVR_2STAGE )
#define LAPACK_EXPORT_ssyevr F77_FUNC( ssyevr , SSYEVR )
#define LAPACK_EXPORT_ssyevx_2stage F77_FUNC( ssyevx_2stage , SSYEVX_2STAGE )
#define LAPACK_EXPORT_ssyevx F77_FUNC( ssyevx , SSYEVX )
#define LAPACK_EXPORT_ssygs2 F77_FUNC( ssygs2 , SSYGS2 )
#define LAPACK_EXPORT_ssygst F77_FUNC( ssygst , SSYGST )
#define LAPACK_EXPORT_ssygv_2stage F77_FUNC( ssygv_2stage , SSYGV_2STAGE )
#define LAPACK_EXPORT_ssygvd F77_FUNC( ssygvd , SSYGVD )
#define LAPACK_EXPORT_ssygv F77_FUNC( ssygv , SSYGV )
#define LAPACK_EXPORT_ssygvx F77_FUNC( ssygvx , SSYGVX )
#define LAPACK_EXPORT_ssyrfs F77_FUNC( ssyrfs , SSYRFS )
#define LAPACK_EXPORT_ssyrfsx F77_FUNC( ssyrfsx , SSYRFSX )
#define LAPACK_EXPORT_ssysv_aa_2stage F77_FUNC( ssysv_aa_2stage , SSYSV_AA_2STAGE )
#define LAPACK_EXPORT_ssysv_aa F77_FUNC( ssysv_aa , SSYSV_AA )
#define LAPACK_EXPORT_ssysv F77_FUNC( ssysv , SSYSV )
#define LAPACK_EXPORT_ssysv_rk F77_FUNC( ssysv_rk , SSYSV_RK )
#define LAPACK_EXPORT_ssysv_rook F77_FUNC( ssysv_rook , SSYSV_ROOK )
#define LAPACK_EXPORT_ssysvx F77_FUNC( ssysvx , SSYSVX )
#define LAPACK_EXPORT_ssysvxx F77_FUNC( ssysvxx , SSYSVXX )
#define LAPACK_EXPORT_ssyswapr F77_FUNC( ssyswapr , SSYSWAPR )
#define LAPACK_EXPORT_ssytd2 F77_FUNC( ssytd2 , SSYTD2 )
#define LAPACK_EXPORT_ssytf2 F77_FUNC( ssytf2 , SSYTF2 )
#define LAPACK_EXPORT_ssytf2_rk F77_FUNC( ssytf2_rk , SSYTF2_RK )
#define LAPACK_EXPORT_ssytf2_rook F77_FUNC( ssytf2_rook , SSYTF2_ROOK )
#define LAPACK_EXPORT_ssytrd_2stage F77_FUNC( ssytrd_2stage , SSYTRD_2STAGE )
#define LAPACK_EXPORT_ssytrd F77_FUNC( ssytrd , SSYTRD )
#define LAPACK_EXPORT_ssytrd_sb2st F77_FUNC( ssytrd_sb2st , SSYTRD_SB2ST )
#define LAPACK_EXPORT_ssytrd_sy2sb F77_FUNC( ssytrd_sy2sb , SSYTRD_SY2SB )
#define LAPACK_EXPORT_ssytrf_aa_2stage F77_FUNC( ssytrf_aa_2stage , SSYTRF_AA_2STAGE )
#define LAPACK_EXPORT_ssytrf_aa F77_FUNC( ssytrf_aa , SSYTRF_AA )
#define LAPACK_EXPORT_ssytrf F77_FUNC( ssytrf , SSYTRF )
#define LAPACK_EXPORT_ssytrf_rk F77_FUNC( ssytrf_rk , SSYTRF_RK )
#define LAPACK_EXPORT_ssytrf_rook F77_FUNC( ssytrf_rook , SSYTRF_ROOK )
#define LAPACK_EXPORT_ssytri2 F77_FUNC( ssytri2 , SSYTRI2 )
#define LAPACK_EXPORT_ssytri2x F77_FUNC( ssytri2x , SSYTRI2X )
#define LAPACK_EXPORT_ssytri_3 F77_FUNC( ssytri_3 , SSYTRI_3 )
#define LAPACK_EXPORT_ssytri_3x F77_FUNC( ssytri_3x , SSYTRI_3X )
#define LAPACK_EXPORT_ssytri F77_FUNC( ssytri , SSYTRI )
#define LAPACK_EXPORT_ssytri_rook F77_FUNC( ssytri_rook , SSYTRI_ROOK )
#define LAPACK_EXPORT_ssytrs2 F77_FUNC( ssytrs2 , SSYTRS2 )
#define LAPACK_EXPORT_ssytrs_3 F77_FUNC( ssytrs_3 , SSYTRS_3 )
#define LAPACK_EXPORT_ssytrs_aa_2stage F77_FUNC( ssytrs_aa_2stage , SSYTRS_AA_2STAGE )
#define LAPACK_EXPORT_ssytrs_aa F77_FUNC( ssytrs_aa , SSYTRS_AA )
#define LAPACK_EXPORT_ssytrs F77_FUNC( ssytrs , SSYTRS )
#define LAPACK_EXPORT_ssytrs_rook F77_FUNC( ssytrs_rook , SSYTRS_ROOK )
#define LAPACK_EXPORT_stbcon F77_FUNC( stbcon , STBCON )
#define LAPACK_EXPORT_stbrfs F77_FUNC( stbrfs , STBRFS )
#define LAPACK_EXPORT_stbtrs F77_FUNC( stbtrs , STBTRS )
#define LAPACK_EXPORT_stfsm F77_FUNC( stfsm , STFSM )
#define LAPACK_EXPORT_stftri F77_FUNC( stftri , STFTRI )
#define LAPACK_EXPORT_stfttp F77_FUNC( stfttp , STFTTP )
#define LAPACK_EXPORT_stfttr F77_FUNC( stfttr , STFTTR )
#define LAPACK_EXPORT_stgevc F77_FUNC( stgevc , STGEVC )
#define LAPACK_EXPORT_stgex2 F77_FUNC( stgex2 , STGEX2 )
#define LAPACK_EXPORT_stgexc F77_FUNC( stgexc , STGEXC )
#define LAPACK_EXPORT_stgsen F77_FUNC( stgsen , STGSEN )
#define LAPACK_EXPORT_stgsja F77_FUNC( stgsja , STGSJA )
#define LAPACK_EXPORT_stgsna F77_FUNC( stgsna , STGSNA )
#define LAPACK_EXPORT_stgsy2 F77_FUNC( stgsy2 , STGSY2 )
#define LAPACK_EXPORT_stgsyl F77_FUNC( stgsyl , STGSYL )
#define LAPACK_EXPORT_stpcon F77_FUNC( stpcon , STPCON )
#define LAPACK_EXPORT_stplqt2 F77_FUNC( stplqt2 , STPLQT2 )
#define LAPACK_EXPORT_stplqt F77_FUNC( stplqt , STPLQT )
#define LAPACK_EXPORT_stpmlqt F77_FUNC( stpmlqt , STPMLQT )
#define LAPACK_EXPORT_stpmqrt F77_FUNC( stpmqrt , STPMQRT )
#define LAPACK_EXPORT_stpqrt2 F77_FUNC( stpqrt2 , STPQRT2 )
#define LAPACK_EXPORT_stpqrt F77_FUNC( stpqrt , STPQRT )
#define LAPACK_EXPORT_stprfb F77_FUNC( stprfb , STPRFB )
#define LAPACK_EXPORT_stprfs F77_FUNC( stprfs , STPRFS )
#define LAPACK_EXPORT_stptri F77_FUNC( stptri , STPTRI )
#define LAPACK_EXPORT_stptrs F77_FUNC( stptrs , STPTRS )
#define LAPACK_EXPORT_stpttf F77_FUNC( stpttf , STPTTF )
#define LAPACK_EXPORT_stpttr F77_FUNC( stpttr , STPTTR )
#define LAPACK_EXPORT_strcon F77_FUNC( strcon , STRCON )
#define LAPACK_EXPORT_strevc3 F77_FUNC( strevc3 , STREVC3 )
#define LAPACK_EXPORT_strevc F77_FUNC( strevc , STREVC )
#define LAPACK_EXPORT_strexc F77_FUNC( strexc , STREXC )
#define LAPACK_EXPORT_strrfs F77_FUNC( strrfs , STRRFS )
#define LAPACK_EXPORT_strsen F77_FUNC( strsen , STRSEN )
#define LAPACK_EXPORT_strsna F77_FUNC( strsna , STRSNA )
#define LAPACK_EXPORT_strsyl F77_FUNC( strsyl , STRSYL )
#define LAPACK_EXPORT_strti2 F77_FUNC( strti2 , STRTI2 )
#define LAPACK_EXPORT_strtri F77_FUNC( strtri , STRTRI )
#define LAPACK_EXPORT_strtrs F77_FUNC( strtrs , STRTRS )
#define LAPACK_EXPORT_strttf F77_FUNC( strttf , STRTTF )
#define LAPACK_EXPORT_strttp F77_FUNC( strttp , STRTTP )
#define LAPACK_EXPORT_stzrqf F77_FUNC( stzrqf , STZRQF )
#define LAPACK_EXPORT_stzrzf F77_FUNC( stzrzf , STZRZF )
#define LAPACK_EXPORT_zbbcsd F77_FUNC( zbbcsd , ZBBCSD )
#define LAPACK_EXPORT_zbdsqr F77_FUNC( zbdsqr , ZBDSQR )
#define LAPACK_EXPORT_zcgesv F77_FUNC( zcgesv , ZCGESV )
#define LAPACK_EXPORT_zcposv F77_FUNC( zcposv , ZCPOSV )
#define LAPACK_EXPORT_zdrscl F77_FUNC( zdrscl , ZDRSCL )
#define LAPACK_EXPORT_zgbbrd F77_FUNC( zgbbrd , ZGBBRD )
#define LAPACK_EXPORT_zgbcon F77_FUNC( zgbcon , ZGBCON )
#define LAPACK_EXPORT_zgbequb F77_FUNC( zgbequb , ZGBEQUB )
#define LAPACK_EXPORT_zgbequ F77_FUNC( zgbequ , ZGBEQU )
#define LAPACK_EXPORT_zgbrfs F77_FUNC( zgbrfs , ZGBRFS )
#define LAPACK_EXPORT_zgbrfsx F77_FUNC( zgbrfsx , ZGBRFSX )
#define LAPACK_EXPORT_zgbsv F77_FUNC( zgbsv , ZGBSV )
#define LAPACK_EXPORT_zgbsvx F77_FUNC( zgbsvx , ZGBSVX )
#define LAPACK_EXPORT_zgbsvxx F77_FUNC( zgbsvxx , ZGBSVXX )
#define LAPACK_EXPORT_zgbtf2 F77_FUNC( zgbtf2 , ZGBTF2 )
#define LAPACK_EXPORT_zgbtrf F77_FUNC( zgbtrf , ZGBTRF )
#define LAPACK_EXPORT_zgbtrs F77_FUNC( zgbtrs , ZGBTRS )
#define LAPACK_EXPORT_zgebak F77_FUNC( zgebak , ZGEBAK )
#define LAPACK_EXPORT_zgebal F77_FUNC( zgebal , ZGEBAL )
#define LAPACK_EXPORT_zgebd2 F77_FUNC( zgebd2 , ZGEBD2 )
#define LAPACK_EXPORT_zgebrd F77_FUNC( zgebrd , ZGEBRD )
#define LAPACK_EXPORT_zgecon F77_FUNC( zgecon , ZGECON )
#define LAPACK_EXPORT_zgeequb F77_FUNC( zgeequb , ZGEEQUB )
#define LAPACK_EXPORT_zgeequ F77_FUNC( zgeequ , ZGEEQU )
#define LAPACK_EXPORT_zgees F77_FUNC( zgees , ZGEES )
#define LAPACK_EXPORT_zgeesx F77_FUNC( zgeesx , ZGEESX )
#define LAPACK_EXPORT_zgeev F77_FUNC( zgeev , ZGEEV )
#define LAPACK_EXPORT_zgeevx F77_FUNC( zgeevx , ZGEEVX )
#define LAPACK_EXPORT_zgegs F77_FUNC( zgegs , ZGEGS )
#define LAPACK_EXPORT_zgegv F77_FUNC( zgegv , ZGEGV )
#define LAPACK_EXPORT_zgehd2 F77_FUNC( zgehd2 , ZGEHD2 )
#define LAPACK_EXPORT_zgehrd F77_FUNC( zgehrd , ZGEHRD )
#define LAPACK_EXPORT_zgejsv F77_FUNC( zgejsv , ZGEJSV )
#define LAPACK_EXPORT_zgelq2 F77_FUNC( zgelq2 , ZGELQ2 )
#define LAPACK_EXPORT_zgelq F77_FUNC( zgelq , ZGELQ )
#define LAPACK_EXPORT_zgelqf F77_FUNC( zgelqf , ZGELQF )
#define LAPACK_EXPORT_zgelqt3 F77_FUNC( zgelqt3 , ZGELQT3 )
#define LAPACK_EXPORT_zgelqt F77_FUNC( zgelqt , ZGELQT )
#define LAPACK_EXPORT_zgelsd F77_FUNC( zgelsd , ZGELSD )
#define LAPACK_EXPORT_zgels F77_FUNC( zgels , ZGELS )
#define LAPACK_EXPORT_zgelss F77_FUNC( zgelss , ZGELSS )
#define LAPACK_EXPORT_zgelsx F77_FUNC( zgelsx , ZGELSX )
#define LAPACK_EXPORT_zgelsy F77_FUNC( zgelsy , ZGELSY )
#define LAPACK_EXPORT_zgemlq F77_FUNC( zgemlq , ZGEMLQ )
#define LAPACK_EXPORT_zgemlqt F77_FUNC( zgemlqt , ZGEMLQT )
#define LAPACK_EXPORT_zgemqr F77_FUNC( zgemqr , ZGEMQR )
#define LAPACK_EXPORT_zgemqrt F77_FUNC( zgemqrt , ZGEMQRT )
#define LAPACK_EXPORT_zgeql2 F77_FUNC( zgeql2 , ZGEQL2 )
#define LAPACK_EXPORT_zgeqlf F77_FUNC( zgeqlf , ZGEQLF )
#define LAPACK_EXPORT_zgeqp3 F77_FUNC( zgeqp3 , ZGEQP3 )
#define LAPACK_EXPORT_zgeqpf F77_FUNC( zgeqpf , ZGEQPF )
#define LAPACK_EXPORT_zgeqr2 F77_FUNC( zgeqr2 , ZGEQR2 )
#define LAPACK_EXPORT_zgeqr2p F77_FUNC( zgeqr2p , ZGEQR2P )
#define LAPACK_EXPORT_zgeqr F77_FUNC( zgeqr , ZGEQR )
#define LAPACK_EXPORT_zgeqrf F77_FUNC( zgeqrf , ZGEQRF )
#define LAPACK_EXPORT_zgeqrfp F77_FUNC( zgeqrfp , ZGEQRFP )
#define LAPACK_EXPORT_zgeqrt2 F77_FUNC( zgeqrt2 , ZGEQRT2 )
#define LAPACK_EXPORT_zgeqrt3 F77_FUNC( zgeqrt3 , ZGEQRT3 )
#define LAPACK_EXPORT_zgeqrt F77_FUNC( zgeqrt , ZGEQRT )
#define LAPACK_EXPORT_zgerfs F77_FUNC( zgerfs , ZGERFS )
#define LAPACK_EXPORT_zgerfsx F77_FUNC( zgerfsx , ZGERFSX )
#define LAPACK_EXPORT_zgerq2 F77_FUNC( zgerq2 , ZGERQ2 )
#define LAPACK_EXPORT_zgerqf F77_FUNC( zgerqf , ZGERQF )
#define LAPACK_EXPORT_zgesc2 F77_FUNC( zgesc2 , ZGESC2 )
#define LAPACK_EXPORT_zgesdd F77_FUNC( zgesdd , ZGESDD )
#define LAPACK_EXPORT_zgesvd F77_FUNC( zgesvd , ZGESVD )
#define LAPACK_EXPORT_zgesvdq F77_FUNC( zgesvdq , ZGESVDQ )
#define LAPACK_EXPORT_zgesvdx F77_FUNC( zgesvdx , ZGESVDX )
#define LAPACK_EXPORT_zgesv F77_FUNC( zgesv , ZGESV )
#define LAPACK_EXPORT_zgesvj F77_FUNC( zgesvj , ZGESVJ )
#define LAPACK_EXPORT_zgesvx F77_FUNC( zgesvx , ZGESVX )
#define LAPACK_EXPORT_zgesvxx F77_FUNC( zgesvxx , ZGESVXX )
#define LAPACK_EXPORT_zgetc2 F77_FUNC( zgetc2 , ZGETC2 )
#define LAPACK_EXPORT_zgetf2 F77_FUNC( zgetf2 , ZGETF2 )
#define LAPACK_EXPORT_zgetrf2 F77_FUNC( zgetrf2 , ZGETRF2 )
#define LAPACK_EXPORT_zgetrf F77_FUNC( zgetrf , ZGETRF )
#define LAPACK_EXPORT_zgetri F77_FUNC( zgetri , ZGETRI )
#define LAPACK_EXPORT_zgetrs F77_FUNC( zgetrs , ZGETRS )
#define LAPACK_EXPORT_zgetsls F77_FUNC( zgetsls , ZGETSLS )
#define LAPACK_EXPORT_zggbak F77_FUNC( zggbak , ZGGBAK )
#define LAPACK_EXPORT_zggbal F77_FUNC( zggbal , ZGGBAL )
#define LAPACK_EXPORT_zgges3 F77_FUNC( zgges3 , ZGGES3 )
#define LAPACK_EXPORT_zgges F77_FUNC( zgges , ZGGES )
#define LAPACK_EXPORT_zggesx F77_FUNC( zggesx , ZGGESX )
#define LAPACK_EXPORT_zggev3 F77_FUNC( zggev3 , ZGGEV3 )
#define LAPACK_EXPORT_zggev F77_FUNC( zggev , ZGGEV )
#define LAPACK_EXPORT_zggevx F77_FUNC( zggevx , ZGGEVX )
#define LAPACK_EXPORT_zggglm F77_FUNC( zggglm , ZGGGLM )
#define LAPACK_EXPORT_zgghd3 F77_FUNC( zgghd3 , ZGGHD3 )
#define LAPACK_EXPORT_zgghrd F77_FUNC( zgghrd , ZGGHRD )
#define LAPACK_EXPORT_zgglse F77_FUNC( zgglse , ZGGLSE )
#define LAPACK_EXPORT_zggqrf F77_FUNC( zggqrf , ZGGQRF )
#define LAPACK_EXPORT_zggrqf F77_FUNC( zggrqf , ZGGRQF )
#define LAPACK_EXPORT_zggsvd3 F77_FUNC( zggsvd3 , ZGGSVD3 )
#define LAPACK_EXPORT_zggsvd F77_FUNC( zggsvd , ZGGSVD )
#define LAPACK_EXPORT_zggsvp3 F77_FUNC( zggsvp3 , ZGGSVP3 )
#define LAPACK_EXPORT_zggsvp F77_FUNC( zggsvp , ZGGSVP )
#define LAPACK_EXPORT_zgsvj0 F77_FUNC( zgsvj0 , ZGSVJ0 )
#define LAPACK_EXPORT_zgsvj1 F77_FUNC( zgsvj1 , ZGSVJ1 )
#define LAPACK_EXPORT_zgtcon F77_FUNC( zgtcon , ZGTCON )
#define LAPACK_EXPORT_zgtrfs F77_FUNC( zgtrfs , ZGTRFS )
#define LAPACK_EXPORT_zgtsv F77_FUNC( zgtsv , ZGTSV )
#define LAPACK_EXPORT_zgtsvx F77_FUNC( zgtsvx , ZGTSVX )
#define LAPACK_EXPORT_zgttrf F77_FUNC( zgttrf , ZGTTRF )
#define LAPACK_EXPORT_zgttrs F77_FUNC( zgttrs , ZGTTRS )
#define LAPACK_EXPORT_zgtts2 F77_FUNC( zgtts2 , ZGTTS2 )
#define LAPACK_EXPORT_zhb2st_kernels F77_FUNC( zhb2st_kernels , ZHB2ST_KERNELS )
#define LAPACK_EXPORT_zhbev_2stage F77_FUNC( zhbev_2stage , ZHBEV_2STAGE )
#define LAPACK_EXPORT_zhbevd_2stage F77_FUNC( zhbevd_2stage , ZHBEVD_2STAGE )
#define LAPACK_EXPORT_zhbevd F77_FUNC( zhbevd , ZHBEVD )
#define LAPACK_EXPORT_zhbev F77_FUNC( zhbev , ZHBEV )
#define LAPACK_EXPORT_zhbevx_2stage F77_FUNC( zhbevx_2stage , ZHBEVX_2STAGE )
#define LAPACK_EXPORT_zhbevx F77_FUNC( zhbevx , ZHBEVX )
#define LAPACK_EXPORT_zhbgst F77_FUNC( zhbgst , ZHBGST )
#define LAPACK_EXPORT_zhbgvd F77_FUNC( zhbgvd , ZHBGVD )
#define LAPACK_EXPORT_zhbgv F77_FUNC( zhbgv , ZHBGV )
#define LAPACK_EXPORT_zhbgvx F77_FUNC( zhbgvx , ZHBGVX )
#define LAPACK_EXPORT_zhbtrd F77_FUNC( zhbtrd , ZHBTRD )
#define LAPACK_EXPORT_zhecon_3 F77_FUNC( zhecon_3 , ZHECON_3 )
#define LAPACK_EXPORT_zhecon F77_FUNC( zhecon , ZHECON )
#define LAPACK_EXPORT_zhecon_rook F77_FUNC( zhecon_rook , ZHECON_ROOK )
#define LAPACK_EXPORT_zheequb F77_FUNC( zheequb , ZHEEQUB )
#define LAPACK_EXPORT_zheev_2stage F77_FUNC( zheev_2stage , ZHEEV_2STAGE )
#define LAPACK_EXPORT_zheevd_2stage F77_FUNC( zheevd_2stage , ZHEEVD_2STAGE )
#define LAPACK_EXPORT_zheevd F77_FUNC( zheevd , ZHEEVD )
#define LAPACK_EXPORT_zheev  F77_FUNC( zheev  , ZHEEV  )
#define LAPACK_EXPORT_zheev F77_FUNC( zheev , ZHEEV )
#define LAPACK_EXPORT_zheevr_2stage F77_FUNC( zheevr_2stage , ZHEEVR_2STAGE )
#define LAPACK_EXPORT_zheevr F77_FUNC( zheevr , ZHEEVR )
#define LAPACK_EXPORT_zheevx_2stage F77_FUNC( zheevx_2stage , ZHEEVX_2STAGE )
#define LAPACK_EXPORT_zheevx F77_FUNC( zheevx , ZHEEVX )
#define LAPACK_EXPORT_zhegs2 F77_FUNC( zhegs2 , ZHEGS2 )
#define LAPACK_EXPORT_zhegst F77_FUNC( zhegst , ZHEGST )
#define LAPACK_EXPORT_zhegv_2stage F77_FUNC( zhegv_2stage , ZHEGV_2STAGE )
#define LAPACK_EXPORT_zhegvd F77_FUNC( zhegvd , ZHEGVD )
#define LAPACK_EXPORT_zhegv F77_FUNC( zhegv , ZHEGV )
#define LAPACK_EXPORT_zhegvx F77_FUNC( zhegvx , ZHEGVX )
#define LAPACK_EXPORT_zherfs F77_FUNC( zherfs , ZHERFS )
#define LAPACK_EXPORT_zherfsx F77_FUNC( zherfsx , ZHERFSX )
#define LAPACK_EXPORT_zhesv_aa_2stage F77_FUNC( zhesv_aa_2stage , ZHESV_AA_2STAGE )
#define LAPACK_EXPORT_zhesv_aa F77_FUNC( zhesv_aa , ZHESV_AA )
#define LAPACK_EXPORT_zhesv F77_FUNC( zhesv , ZHESV )
#define LAPACK_EXPORT_zhesv_rk F77_FUNC( zhesv_rk , ZHESV_RK )
#define LAPACK_EXPORT_zhesv_rook F77_FUNC( zhesv_rook , ZHESV_ROOK )
#define LAPACK_EXPORT_zhesvx F77_FUNC( zhesvx , ZHESVX )
#define LAPACK_EXPORT_zhesvxx F77_FUNC( zhesvxx , ZHESVXX )
#define LAPACK_EXPORT_zheswapr F77_FUNC( zheswapr , ZHESWAPR )
#define LAPACK_EXPORT_zhetd2 F77_FUNC( zhetd2 , ZHETD2 )
#define LAPACK_EXPORT_zhetf2 F77_FUNC( zhetf2 , ZHETF2 )
#define LAPACK_EXPORT_zhetf2_rk F77_FUNC( zhetf2_rk , ZHETF2_RK )
#define LAPACK_EXPORT_zhetf2_rook F77_FUNC( zhetf2_rook , ZHETF2_ROOK )
#define LAPACK_EXPORT_zhetrd_2stage F77_FUNC( zhetrd_2stage , ZHETRD_2STAGE )
#define LAPACK_EXPORT_zhetrd F77_FUNC( zhetrd , ZHETRD )
#define LAPACK_EXPORT_zhetrd_hb2st F77_FUNC( zhetrd_hb2st , ZHETRD_HB2ST )
#define LAPACK_EXPORT_zhetrd_he2hb F77_FUNC( zhetrd_he2hb , ZHETRD_HE2HB )
#define LAPACK_EXPORT_zhetrf_aa_2stage F77_FUNC( zhetrf_aa_2stage , ZHETRF_AA_2STAGE )
#define LAPACK_EXPORT_zhetrf_aa F77_FUNC( zhetrf_aa , ZHETRF_AA )
#define LAPACK_EXPORT_zhetrf F77_FUNC( zhetrf , ZHETRF )
#define LAPACK_EXPORT_zhetrf_rk F77_FUNC( zhetrf_rk , ZHETRF_RK )
#define LAPACK_EXPORT_zhetrf_rook F77_FUNC( zhetrf_rook , ZHETRF_ROOK )
#define LAPACK_EXPORT_zhetri2 F77_FUNC( zhetri2 , ZHETRI2 )
#define LAPACK_EXPORT_zhetri2x F77_FUNC( zhetri2x , ZHETRI2X )
#define LAPACK_EXPORT_zhetri_3 F77_FUNC( zhetri_3 , ZHETRI_3 )
#define LAPACK_EXPORT_zhetri_3x F77_FUNC( zhetri_3x , ZHETRI_3X )
#define LAPACK_EXPORT_zhetri F77_FUNC( zhetri , ZHETRI )
#define LAPACK_EXPORT_zhetri_rook F77_FUNC( zhetri_rook , ZHETRI_ROOK )
#define LAPACK_EXPORT_zhetrs2 F77_FUNC( zhetrs2 , ZHETRS2 )
#define LAPACK_EXPORT_zhetrs_3 F77_FUNC( zhetrs_3 , ZHETRS_3 )
#define LAPACK_EXPORT_zhetrs_aa_2stage F77_FUNC( zhetrs_aa_2stage , ZHETRS_AA_2STAGE )
#define LAPACK_EXPORT_zhetrs_aa F77_FUNC( zhetrs_aa , ZHETRS_AA )
#define LAPACK_EXPORT_zhetrs F77_FUNC( zhetrs , ZHETRS )
#define LAPACK_EXPORT_zhetrs_rook F77_FUNC( zhetrs_rook , ZHETRS_ROOK )
#define LAPACK_EXPORT_zhfrk F77_FUNC( zhfrk , ZHFRK )
#define LAPACK_EXPORT_zhgeqz F77_FUNC( zhgeqz , ZHGEQZ )
#define LAPACK_EXPORT_zhpcon F77_FUNC( zhpcon , ZHPCON )
#define LAPACK_EXPORT_zhpevd F77_FUNC( zhpevd , ZHPEVD )
#define LAPACK_EXPORT_zhpev F77_FUNC( zhpev , ZHPEV )
#define LAPACK_EXPORT_zhpevx F77_FUNC( zhpevx , ZHPEVX )
#define LAPACK_EXPORT_zhpgst F77_FUNC( zhpgst , ZHPGST )
#define LAPACK_EXPORT_zhpgvd F77_FUNC( zhpgvd , ZHPGVD )
#define LAPACK_EXPORT_zhpgv F77_FUNC( zhpgv , ZHPGV )
#define LAPACK_EXPORT_zhpgvx F77_FUNC( zhpgvx , ZHPGVX )
#define LAPACK_EXPORT_zhprfs F77_FUNC( zhprfs , ZHPRFS )
#define LAPACK_EXPORT_zhpsv F77_FUNC( zhpsv , ZHPSV )
#define LAPACK_EXPORT_zhpsvx F77_FUNC( zhpsvx , ZHPSVX )
#define LAPACK_EXPORT_zhptrd F77_FUNC( zhptrd , ZHPTRD )
#define LAPACK_EXPORT_zhptrf F77_FUNC( zhptrf , ZHPTRF )
#define LAPACK_EXPORT_zhptri F77_FUNC( zhptri , ZHPTRI )
#define LAPACK_EXPORT_zhptrs F77_FUNC( zhptrs , ZHPTRS )
#define LAPACK_EXPORT_zhsein F77_FUNC( zhsein , ZHSEIN )
#define LAPACK_EXPORT_zhseqr F77_FUNC( zhseqr , ZHSEQR )
#define LAPACK_EXPORT_zlabrd F77_FUNC( zlabrd , ZLABRD )
#define LAPACK_EXPORT_zlacgv F77_FUNC( zlacgv , ZLACGV )
#define LAPACK_EXPORT_zlacn2 F77_FUNC( zlacn2 , ZLACN2 )
#define LAPACK_EXPORT_zlacon F77_FUNC( zlacon , ZLACON )
#define LAPACK_EXPORT_zlacp2 F77_FUNC( zlacp2 , ZLACP2 )
#define LAPACK_EXPORT_zlacpy F77_FUNC( zlacpy , ZLACPY )
#define LAPACK_EXPORT_zlacrm F77_FUNC( zlacrm , ZLACRM )
#define LAPACK_EXPORT_zlacrt F77_FUNC( zlacrt , ZLACRT )
#define LAPACK_EXPORT_zladiv F77_FUNC( zladiv , ZLADIV )
#define LAPACK_EXPORT_zlaed0 F77_FUNC( zlaed0 , ZLAED0 )
#define LAPACK_EXPORT_zlaed7 F77_FUNC( zlaed7 , ZLAED7 )
#define LAPACK_EXPORT_zlaed8 F77_FUNC( zlaed8 , ZLAED8 )
#define LAPACK_EXPORT_zlaein F77_FUNC( zlaein , ZLAEIN )
#define LAPACK_EXPORT_zlaesy F77_FUNC( zlaesy , ZLAESY )
#define LAPACK_EXPORT_zlaev2 F77_FUNC( zlaev2 , ZLAEV2 )
#define LAPACK_EXPORT_zlag2c F77_FUNC( zlag2c , ZLAG2C )
#define LAPACK_EXPORT_zla_gbamv F77_FUNC( zla_gbamv , ZLA_GBAMV )
#define LAPACK_EXPORT_zla_gbrcond_c F77_FUNC( zla_gbrcond_c , ZLA_GBRCOND_C )
#define LAPACK_EXPORT_zla_gbrcond_x F77_FUNC( zla_gbrcond_x , ZLA_GBRCOND_X )
#define LAPACK_EXPORT_zla_gbrfsx_extended F77_FUNC( zla_gbrfsx_extended , ZLA_GBRFSX_EXTENDED )
#define LAPACK_EXPORT_zla_gbrpvgrw F77_FUNC( zla_gbrpvgrw , ZLA_GBRPVGRW )
#define LAPACK_EXPORT_zla_geamv F77_FUNC( zla_geamv , ZLA_GEAMV )
#define LAPACK_EXPORT_zla_gercond_c F77_FUNC( zla_gercond_c , ZLA_GERCOND_C )
#define LAPACK_EXPORT_zla_gercond_x F77_FUNC( zla_gercond_x , ZLA_GERCOND_X )
#define LAPACK_EXPORT_zla_gerfsx_extended F77_FUNC( zla_gerfsx_extended , ZLA_GERFSX_EXTENDED )
#define LAPACK_EXPORT_zla_gerpvgrw F77_FUNC( zla_gerpvgrw , ZLA_GERPVGRW )
#define LAPACK_EXPORT_zlags2 F77_FUNC( zlags2 , ZLAGS2 )
#define LAPACK_EXPORT_zlagtm F77_FUNC( zlagtm , ZLAGTM )
#define LAPACK_EXPORT_zla_heamv F77_FUNC( zla_heamv , ZLA_HEAMV )
#define LAPACK_EXPORT_zlahef_aa F77_FUNC( zlahef_aa , ZLAHEF_AA )
#define LAPACK_EXPORT_zlahef F77_FUNC( zlahef , ZLAHEF )
#define LAPACK_EXPORT_zlahef_rk F77_FUNC( zlahef_rk , ZLAHEF_RK )
#define LAPACK_EXPORT_zlahef_rook F77_FUNC( zlahef_rook , ZLAHEF_ROOK )
#define LAPACK_EXPORT_zla_hercond_c F77_FUNC( zla_hercond_c , ZLA_HERCOND_C )
#define LAPACK_EXPORT_zla_hercond_x F77_FUNC( zla_hercond_x , ZLA_HERCOND_X )
#define LAPACK_EXPORT_zla_herfsx_extended F77_FUNC( zla_herfsx_extended , ZLA_HERFSX_EXTENDED )
#define LAPACK_EXPORT_zla_herpvgrw F77_FUNC( zla_herpvgrw , ZLA_HERPVGRW )
#define LAPACK_EXPORT_zlahqr F77_FUNC( zlahqr , ZLAHQR )
#define LAPACK_EXPORT_zlahr2 F77_FUNC( zlahr2 , ZLAHR2 )
#define LAPACK_EXPORT_zlahrd F77_FUNC( zlahrd , ZLAHRD )
#define LAPACK_EXPORT_zlaic1 F77_FUNC( zlaic1 , ZLAIC1 )
#define LAPACK_EXPORT_zla_lin_berr F77_FUNC( zla_lin_berr , ZLA_LIN_BERR )
#define LAPACK_EXPORT_zlals0 F77_FUNC( zlals0 , ZLALS0 )
#define LAPACK_EXPORT_zlalsa F77_FUNC( zlalsa , ZLALSA )
#define LAPACK_EXPORT_zlalsd F77_FUNC( zlalsd , ZLALSD )
#define LAPACK_EXPORT_zlamswlq F77_FUNC( zlamswlq , ZLAMSWLQ )
#define LAPACK_EXPORT_zlamtsqr F77_FUNC( zlamtsqr , ZLAMTSQR )
#define LAPACK_EXPORT_zlangb F77_FUNC( zlangb , ZLANGB )
#define LAPACK_EXPORT_zlange F77_FUNC( zlange , ZLANGE )
#define LAPACK_EXPORT_zlangt F77_FUNC( zlangt , ZLANGT )
#define LAPACK_EXPORT_zlanhb F77_FUNC( zlanhb , ZLANHB )
#define LAPACK_EXPORT_zlanhe F77_FUNC( zlanhe , ZLANHE )
#define LAPACK_EXPORT_zlanhf F77_FUNC( zlanhf , ZLANHF )
#define LAPACK_EXPORT_zlanhp F77_FUNC( zlanhp , ZLANHP )
#define LAPACK_EXPORT_zlanhs F77_FUNC( zlanhs , ZLANHS )
#define LAPACK_EXPORT_zlanht F77_FUNC( zlanht , ZLANHT )
#define LAPACK_EXPORT_zlansb F77_FUNC( zlansb , ZLANSB )
#define LAPACK_EXPORT_zlansp F77_FUNC( zlansp , ZLANSP )
#define LAPACK_EXPORT_zlansy F77_FUNC( zlansy , ZLANSY )
#define LAPACK_EXPORT_zlantb F77_FUNC( zlantb , ZLANTB )
#define LAPACK_EXPORT_zlantp F77_FUNC( zlantp , ZLANTP )
#define LAPACK_EXPORT_zlantr F77_FUNC( zlantr , ZLANTR )
#define LAPACK_EXPORT_zlapll F77_FUNC( zlapll , ZLAPLL )
#define LAPACK_EXPORT_zlapmr F77_FUNC( zlapmr , ZLAPMR )
#define LAPACK_EXPORT_zlapmt F77_FUNC( zlapmt , ZLAPMT )
#define LAPACK_EXPORT_zla_porcond_c F77_FUNC( zla_porcond_c , ZLA_PORCOND_C )
#define LAPACK_EXPORT_zla_porcond_x F77_FUNC( zla_porcond_x , ZLA_PORCOND_X )
#define LAPACK_EXPORT_zla_porfsx_extended F77_FUNC( zla_porfsx_extended , ZLA_PORFSX_EXTENDED )
#define LAPACK_EXPORT_zla_porpvgrw F77_FUNC( zla_porpvgrw , ZLA_PORPVGRW )
#define LAPACK_EXPORT_zlaqgb F77_FUNC( zlaqgb , ZLAQGB )
#define LAPACK_EXPORT_zlaqge F77_FUNC( zlaqge , ZLAQGE )
#define LAPACK_EXPORT_zlaqhb F77_FUNC( zlaqhb , ZLAQHB )
#define LAPACK_EXPORT_zlaqhe F77_FUNC( zlaqhe , ZLAQHE )
#define LAPACK_EXPORT_zlaqhp F77_FUNC( zlaqhp , ZLAQHP )
#define LAPACK_EXPORT_zlaqp2 F77_FUNC( zlaqp2 , ZLAQP2 )
#define LAPACK_EXPORT_zlaqps F77_FUNC( zlaqps , ZLAQPS )
#define LAPACK_EXPORT_zlaqr0 F77_FUNC( zlaqr0 , ZLAQR0 )
#define LAPACK_EXPORT_zlaqr1 F77_FUNC( zlaqr1 , ZLAQR1 )
#define LAPACK_EXPORT_zlaqr2 F77_FUNC( zlaqr2 , ZLAQR2 )
#define LAPACK_EXPORT_zlaqr3 F77_FUNC( zlaqr3 , ZLAQR3 )
#define LAPACK_EXPORT_zlaqr4 F77_FUNC( zlaqr4 , ZLAQR4 )
#define LAPACK_EXPORT_zlaqr5 F77_FUNC( zlaqr5 , ZLAQR5 )
#define LAPACK_EXPORT_zlaqsb F77_FUNC( zlaqsb , ZLAQSB )
#define LAPACK_EXPORT_zlaqsp F77_FUNC( zlaqsp , ZLAQSP )
#define LAPACK_EXPORT_zlaqsy F77_FUNC( zlaqsy , ZLAQSY )
#define LAPACK_EXPORT_zlar1v F77_FUNC( zlar1v , ZLAR1V )
#define LAPACK_EXPORT_zlar2v F77_FUNC( zlar2v , ZLAR2V )
#define LAPACK_EXPORT_zlarcm F77_FUNC( zlarcm , ZLARCM )
#define LAPACK_EXPORT_zlarfb F77_FUNC( zlarfb , ZLARFB )
#define LAPACK_EXPORT_zlarf F77_FUNC( zlarf , ZLARF )
#define LAPACK_EXPORT_zlarfg F77_FUNC( zlarfg , ZLARFG )
#define LAPACK_EXPORT_zlarfgp F77_FUNC( zlarfgp , ZLARFGP )
#define LAPACK_EXPORT_zlarft F77_FUNC( zlarft , ZLARFT )
#define LAPACK_EXPORT_zlarfx F77_FUNC( zlarfx , ZLARFX )
#define LAPACK_EXPORT_zlarfy F77_FUNC( zlarfy , ZLARFY )
#define LAPACK_EXPORT_zlargv F77_FUNC( zlargv , ZLARGV )
#define LAPACK_EXPORT_zlarnv F77_FUNC( zlarnv , ZLARNV )
#define LAPACK_EXPORT_zlarrv F77_FUNC( zlarrv , ZLARRV )
#define LAPACK_EXPORT_zlarscl2 F77_FUNC( zlarscl2 , ZLARSCL2 )
#define LAPACK_EXPORT_zlartg F77_FUNC( zlartg , ZLARTG )
#define LAPACK_EXPORT_zlartv F77_FUNC( zlartv , ZLARTV )
#define LAPACK_EXPORT_zlarzb F77_FUNC( zlarzb , ZLARZB )
#define LAPACK_EXPORT_zlarz F77_FUNC( zlarz , ZLARZ )
#define LAPACK_EXPORT_zlarzt F77_FUNC( zlarzt , ZLARZT )
#define LAPACK_EXPORT_zlascl2 F77_FUNC( zlascl2 , ZLASCL2 )
#define LAPACK_EXPORT_zlascl F77_FUNC( zlascl , ZLASCL )
#define LAPACK_EXPORT_zlaset F77_FUNC( zlaset , ZLASET )
#define LAPACK_EXPORT_zlasr F77_FUNC( zlasr , ZLASR )
#define LAPACK_EXPORT_zlassq F77_FUNC( zlassq , ZLASSQ )
#define LAPACK_EXPORT_zlaswlq F77_FUNC( zlaswlq , ZLASWLQ )
#define LAPACK_EXPORT_zlaswp F77_FUNC( zlaswp , ZLASWP )
#define LAPACK_EXPORT_zla_syamv F77_FUNC( zla_syamv , ZLA_SYAMV )
#define LAPACK_EXPORT_zlasyf_aa F77_FUNC( zlasyf_aa , ZLASYF_AA )
#define LAPACK_EXPORT_zlasyf F77_FUNC( zlasyf , ZLASYF )
#define LAPACK_EXPORT_zlasyf_rk F77_FUNC( zlasyf_rk , ZLASYF_RK )
#define LAPACK_EXPORT_zlasyf_rook F77_FUNC( zlasyf_rook , ZLASYF_ROOK )
#define LAPACK_EXPORT_zla_syrcond_c F77_FUNC( zla_syrcond_c , ZLA_SYRCOND_C )
#define LAPACK_EXPORT_zla_syrcond_x F77_FUNC( zla_syrcond_x , ZLA_SYRCOND_X )
#define LAPACK_EXPORT_zla_syrfsx_extended F77_FUNC( zla_syrfsx_extended , ZLA_SYRFSX_EXTENDED )
#define LAPACK_EXPORT_zla_syrpvgrw F77_FUNC( zla_syrpvgrw , ZLA_SYRPVGRW )
#define LAPACK_EXPORT_zlat2c F77_FUNC( zlat2c , ZLAT2C )
#define LAPACK_EXPORT_zlatbs F77_FUNC( zlatbs , ZLATBS )
#define LAPACK_EXPORT_zlatdf F77_FUNC( zlatdf , ZLATDF )
#define LAPACK_EXPORT_zlatps F77_FUNC( zlatps , ZLATPS )
#define LAPACK_EXPORT_zlatrd F77_FUNC( zlatrd , ZLATRD )
#define LAPACK_EXPORT_zlatrs F77_FUNC( zlatrs , ZLATRS )
#define LAPACK_EXPORT_zlatrz F77_FUNC( zlatrz , ZLATRZ )
#define LAPACK_EXPORT_zlatsqr F77_FUNC( zlatsqr , ZLATSQR )
#define LAPACK_EXPORT_zlatzm F77_FUNC( zlatzm , ZLATZM )
#define LAPACK_EXPORT_zlaunhr_col_getrfnp2 F77_FUNC( zlaunhr_col_getrfnp2 , ZLAUNHR_COL_GETRFNP2 )
#define LAPACK_EXPORT_zlaunhr_col_getrfnp F77_FUNC( zlaunhr_col_getrfnp , ZLAUNHR_COL_GETRFNP )
#define LAPACK_EXPORT_zlauu2 F77_FUNC( zlauu2 , ZLAUU2 )
#define LAPACK_EXPORT_zlauum F77_FUNC( zlauum , ZLAUUM )
#define LAPACK_EXPORT_zla_wwaddw F77_FUNC( zla_wwaddw , ZLA_WWADDW )
#define LAPACK_EXPORT_zpbcon F77_FUNC( zpbcon , ZPBCON )
#define LAPACK_EXPORT_zpbequ F77_FUNC( zpbequ , ZPBEQU )
#define LAPACK_EXPORT_zpbrfs F77_FUNC( zpbrfs , ZPBRFS )
#define LAPACK_EXPORT_zpbstf F77_FUNC( zpbstf , ZPBSTF )
#define LAPACK_EXPORT_zpbsv F77_FUNC( zpbsv , ZPBSV )
#define LAPACK_EXPORT_zpbsvx F77_FUNC( zpbsvx , ZPBSVX )
#define LAPACK_EXPORT_zpbtf2 F77_FUNC( zpbtf2 , ZPBTF2 )
#define LAPACK_EXPORT_zpbtrf F77_FUNC( zpbtrf , ZPBTRF )
#define LAPACK_EXPORT_zpbtrs F77_FUNC( zpbtrs , ZPBTRS )
#define LAPACK_EXPORT_zpftrf F77_FUNC( zpftrf , ZPFTRF )
#define LAPACK_EXPORT_zpftri F77_FUNC( zpftri , ZPFTRI )
#define LAPACK_EXPORT_zpftrs F77_FUNC( zpftrs , ZPFTRS )
#define LAPACK_EXPORT_zpocon F77_FUNC( zpocon , ZPOCON )
#define LAPACK_EXPORT_zpoequb F77_FUNC( zpoequb , ZPOEQUB )
#define LAPACK_EXPORT_zpoequ F77_FUNC( zpoequ , ZPOEQU )
#define LAPACK_EXPORT_zporfs F77_FUNC( zporfs , ZPORFS )
#define LAPACK_EXPORT_zporfsx F77_FUNC( zporfsx , ZPORFSX )
#define LAPACK_EXPORT_zposv F77_FUNC( zposv , ZPOSV )
#define LAPACK_EXPORT_zposvx F77_FUNC( zposvx , ZPOSVX )
#define LAPACK_EXPORT_zposvxx F77_FUNC( zposvxx , ZPOSVXX )
#define LAPACK_EXPORT_zpotf2 F77_FUNC( zpotf2 , ZPOTF2 )
#define LAPACK_EXPORT_zpotrf2 F77_FUNC( zpotrf2 , ZPOTRF2 )
#define LAPACK_EXPORT_zpotrf F77_FUNC( zpotrf , ZPOTRF )
#define LAPACK_EXPORT_zpotri F77_FUNC( zpotri , ZPOTRI )
#define LAPACK_EXPORT_zpotrs F77_FUNC( zpotrs , ZPOTRS )
#define LAPACK_EXPORT_zppcon F77_FUNC( zppcon , ZPPCON )
#define LAPACK_EXPORT_zppequ F77_FUNC( zppequ , ZPPEQU )
#define LAPACK_EXPORT_zpprfs F77_FUNC( zpprfs , ZPPRFS )
#define LAPACK_EXPORT_zppsv F77_FUNC( zppsv , ZPPSV )
#define LAPACK_EXPORT_zppsvx F77_FUNC( zppsvx , ZPPSVX )
#define LAPACK_EXPORT_zpptrf F77_FUNC( zpptrf , ZPPTRF )
#define LAPACK_EXPORT_zpptri F77_FUNC( zpptri , ZPPTRI )
#define LAPACK_EXPORT_zpptrs F77_FUNC( zpptrs , ZPPTRS )
#define LAPACK_EXPORT_zpstf2 F77_FUNC( zpstf2 , ZPSTF2 )
#define LAPACK_EXPORT_zpstrf F77_FUNC( zpstrf , ZPSTRF )
#define LAPACK_EXPORT_zptcon F77_FUNC( zptcon , ZPTCON )
#define LAPACK_EXPORT_zpteqr F77_FUNC( zpteqr , ZPTEQR )
#define LAPACK_EXPORT_zptrfs F77_FUNC( zptrfs , ZPTRFS )
#define LAPACK_EXPORT_zptsv F77_FUNC( zptsv , ZPTSV )
#define LAPACK_EXPORT_zptsvx F77_FUNC( zptsvx , ZPTSVX )
#define LAPACK_EXPORT_zpttrf F77_FUNC( zpttrf , ZPTTRF )
#define LAPACK_EXPORT_zpttrs F77_FUNC( zpttrs , ZPTTRS )
#define LAPACK_EXPORT_zptts2 F77_FUNC( zptts2 , ZPTTS2 )
#define LAPACK_EXPORT_zrot F77_FUNC( zrot , ZROT )
#define LAPACK_EXPORT_zspcon F77_FUNC( zspcon , ZSPCON )
#define LAPACK_EXPORT_zspmv F77_FUNC( zspmv , ZSPMV )
#define LAPACK_EXPORT_zspr F77_FUNC( zspr , ZSPR )
#define LAPACK_EXPORT_zsprfs F77_FUNC( zsprfs , ZSPRFS )
#define LAPACK_EXPORT_zspsv F77_FUNC( zspsv , ZSPSV )
#define LAPACK_EXPORT_zspsvx F77_FUNC( zspsvx , ZSPSVX )
#define LAPACK_EXPORT_zsptrf F77_FUNC( zsptrf , ZSPTRF )
#define LAPACK_EXPORT_zsptri F77_FUNC( zsptri , ZSPTRI )
#define LAPACK_EXPORT_zsptrs F77_FUNC( zsptrs , ZSPTRS )
#define LAPACK_EXPORT_zstedc F77_FUNC( zstedc , ZSTEDC )
#define LAPACK_EXPORT_zstegr F77_FUNC( zstegr , ZSTEGR )
#define LAPACK_EXPORT_zstein F77_FUNC( zstein , ZSTEIN )
#define LAPACK_EXPORT_zstemr F77_FUNC( zstemr , ZSTEMR )
#define LAPACK_EXPORT_zsteqr F77_FUNC( zsteqr , ZSTEQR )
#define LAPACK_EXPORT_zsycon_3 F77_FUNC( zsycon_3 , ZSYCON_3 )
#define LAPACK_EXPORT_zsycon F77_FUNC( zsycon , ZSYCON )
#define LAPACK_EXPORT_zsycon_rook F77_FUNC( zsycon_rook , ZSYCON_ROOK )
#define LAPACK_EXPORT_zsyconv F77_FUNC( zsyconv , ZSYCONV )
#define LAPACK_EXPORT_zsyconvf F77_FUNC( zsyconvf , ZSYCONVF )
#define LAPACK_EXPORT_zsyconvf_rook F77_FUNC( zsyconvf_rook , ZSYCONVF_ROOK )
#define LAPACK_EXPORT_zsyequb F77_FUNC( zsyequb , ZSYEQUB )
#define LAPACK_EXPORT_zsymv F77_FUNC( zsymv , ZSYMV )
#define LAPACK_EXPORT_zsyr F77_FUNC( zsyr , ZSYR )
#define LAPACK_EXPORT_zsyrfs F77_FUNC( zsyrfs , ZSYRFS )
#define LAPACK_EXPORT_zsyrfsx F77_FUNC( zsyrfsx , ZSYRFSX )
#define LAPACK_EXPORT_zsysv_aa_2stage F77_FUNC( zsysv_aa_2stage , ZSYSV_AA_2STAGE )
#define LAPACK_EXPORT_zsysv_aa F77_FUNC( zsysv_aa , ZSYSV_AA )
#define LAPACK_EXPORT_zsysv F77_FUNC( zsysv , ZSYSV )
#define LAPACK_EXPORT_zsysv_rk F77_FUNC( zsysv_rk , ZSYSV_RK )
#define LAPACK_EXPORT_zsysv_rook F77_FUNC( zsysv_rook , ZSYSV_ROOK )
#define LAPACK_EXPORT_zsysvx F77_FUNC( zsysvx , ZSYSVX )
#define LAPACK_EXPORT_zsysvxx F77_FUNC( zsysvxx , ZSYSVXX )
#define LAPACK_EXPORT_zsyswapr F77_FUNC( zsyswapr , ZSYSWAPR )
#define LAPACK_EXPORT_zsytf2 F77_FUNC( zsytf2 , ZSYTF2 )
#define LAPACK_EXPORT_zsytf2_rk F77_FUNC( zsytf2_rk , ZSYTF2_RK )
#define LAPACK_EXPORT_zsytf2_rook F77_FUNC( zsytf2_rook , ZSYTF2_ROOK )
#define LAPACK_EXPORT_zsytrf_aa_2stage F77_FUNC( zsytrf_aa_2stage , ZSYTRF_AA_2STAGE )
#define LAPACK_EXPORT_zsytrf_aa F77_FUNC( zsytrf_aa , ZSYTRF_AA )
#define LAPACK_EXPORT_zsytrf F77_FUNC( zsytrf , ZSYTRF )
#define LAPACK_EXPORT_zsytrf_rk F77_FUNC( zsytrf_rk , ZSYTRF_RK )
#define LAPACK_EXPORT_zsytrf_rook F77_FUNC( zsytrf_rook , ZSYTRF_ROOK )
#define LAPACK_EXPORT_zsytri2 F77_FUNC( zsytri2 , ZSYTRI2 )
#define LAPACK_EXPORT_zsytri2x F77_FUNC( zsytri2x , ZSYTRI2X )
#define LAPACK_EXPORT_zsytri_3 F77_FUNC( zsytri_3 , ZSYTRI_3 )
#define LAPACK_EXPORT_zsytri_3x F77_FUNC( zsytri_3x , ZSYTRI_3X )
#define LAPACK_EXPORT_zsytri F77_FUNC( zsytri , ZSYTRI )
#define LAPACK_EXPORT_zsytri_rook F77_FUNC( zsytri_rook , ZSYTRI_ROOK )
#define LAPACK_EXPORT_zsytrs2 F77_FUNC( zsytrs2 , ZSYTRS2 )
#define LAPACK_EXPORT_zsytrs_3 F77_FUNC( zsytrs_3 , ZSYTRS_3 )
#define LAPACK_EXPORT_zsytrs_aa_2stage F77_FUNC( zsytrs_aa_2stage , ZSYTRS_AA_2STAGE )
#define LAPACK_EXPORT_zsytrs_aa F77_FUNC( zsytrs_aa , ZSYTRS_AA )
#define LAPACK_EXPORT_zsytrs F77_FUNC( zsytrs , ZSYTRS )
#define LAPACK_EXPORT_zsytrs_rook F77_FUNC( zsytrs_rook , ZSYTRS_ROOK )
#define LAPACK_EXPORT_ztbcon F77_FUNC( ztbcon , ZTBCON )
#define LAPACK_EXPORT_ztbrfs F77_FUNC( ztbrfs , ZTBRFS )
#define LAPACK_EXPORT_ztbtrs F77_FUNC( ztbtrs , ZTBTRS )
#define LAPACK_EXPORT_ztfsm F77_FUNC( ztfsm , ZTFSM )
#define LAPACK_EXPORT_ztftri F77_FUNC( ztftri , ZTFTRI )
#define LAPACK_EXPORT_ztfttp F77_FUNC( ztfttp , ZTFTTP )
#define LAPACK_EXPORT_ztfttr F77_FUNC( ztfttr , ZTFTTR )
#define LAPACK_EXPORT_ztgevc F77_FUNC( ztgevc , ZTGEVC )
#define LAPACK_EXPORT_ztgex2 F77_FUNC( ztgex2 , ZTGEX2 )
#define LAPACK_EXPORT_ztgexc F77_FUNC( ztgexc , ZTGEXC )
#define LAPACK_EXPORT_ztgsen F77_FUNC( ztgsen , ZTGSEN )
#define LAPACK_EXPORT_ztgsja F77_FUNC( ztgsja , ZTGSJA )
#define LAPACK_EXPORT_ztgsna F77_FUNC( ztgsna , ZTGSNA )
#define LAPACK_EXPORT_ztgsy2 F77_FUNC( ztgsy2 , ZTGSY2 )
#define LAPACK_EXPORT_ztgsyl F77_FUNC( ztgsyl , ZTGSYL )
#define LAPACK_EXPORT_ztpcon F77_FUNC( ztpcon , ZTPCON )
#define LAPACK_EXPORT_ztplqt2 F77_FUNC( ztplqt2 , ZTPLQT2 )
#define LAPACK_EXPORT_ztplqt F77_FUNC( ztplqt , ZTPLQT )
#define LAPACK_EXPORT_ztpmlqt F77_FUNC( ztpmlqt , ZTPMLQT )
#define LAPACK_EXPORT_ztpmqrt F77_FUNC( ztpmqrt , ZTPMQRT )
#define LAPACK_EXPORT_ztpqrt2 F77_FUNC( ztpqrt2 , ZTPQRT2 )
#define LAPACK_EXPORT_ztpqrt F77_FUNC( ztpqrt , ZTPQRT )
#define LAPACK_EXPORT_ztprfb F77_FUNC( ztprfb , ZTPRFB )
#define LAPACK_EXPORT_ztprfs F77_FUNC( ztprfs , ZTPRFS )
#define LAPACK_EXPORT_ztptri F77_FUNC( ztptri , ZTPTRI )
#define LAPACK_EXPORT_ztptrs F77_FUNC( ztptrs , ZTPTRS )
#define LAPACK_EXPORT_ztpttf F77_FUNC( ztpttf , ZTPTTF )
#define LAPACK_EXPORT_ztpttr F77_FUNC( ztpttr , ZTPTTR )
#define LAPACK_EXPORT_ztrcon F77_FUNC( ztrcon , ZTRCON )
#define LAPACK_EXPORT_ztrevc3 F77_FUNC( ztrevc3 , ZTREVC3 )
#define LAPACK_EXPORT_ztrevc F77_FUNC( ztrevc , ZTREVC )
#define LAPACK_EXPORT_ztrexc F77_FUNC( ztrexc , ZTREXC )
#define LAPACK_EXPORT_ztrrfs F77_FUNC( ztrrfs , ZTRRFS )
#define LAPACK_EXPORT_ztrsen F77_FUNC( ztrsen , ZTRSEN )
#define LAPACK_EXPORT_ztrsna F77_FUNC( ztrsna , ZTRSNA )
#define LAPACK_EXPORT_ztrsyl F77_FUNC( ztrsyl , ZTRSYL )
#define LAPACK_EXPORT_ztrti2 F77_FUNC( ztrti2 , ZTRTI2 )
#define LAPACK_EXPORT_ztrtri F77_FUNC( ztrtri , ZTRTRI )
#define LAPACK_EXPORT_ztrtrs F77_FUNC( ztrtrs , ZTRTRS )
#define LAPACK_EXPORT_ztrttf F77_FUNC( ztrttf , ZTRTTF )
#define LAPACK_EXPORT_ztrttp F77_FUNC( ztrttp , ZTRTTP )
#define LAPACK_EXPORT_ztzrqf F77_FUNC( ztzrqf , ZTZRQF )
#define LAPACK_EXPORT_ztzrzf F77_FUNC( ztzrzf , ZTZRZF )
#define LAPACK_EXPORT_zunbdb1 F77_FUNC( zunbdb1 , ZUNBDB1 )
#define LAPACK_EXPORT_zunbdb2 F77_FUNC( zunbdb2 , ZUNBDB2 )
#define LAPACK_EXPORT_zunbdb3 F77_FUNC( zunbdb3 , ZUNBDB3 )
#define LAPACK_EXPORT_zunbdb4 F77_FUNC( zunbdb4 , ZUNBDB4 )
#define LAPACK_EXPORT_zunbdb5 F77_FUNC( zunbdb5 , ZUNBDB5 )
#define LAPACK_EXPORT_zunbdb6 F77_FUNC( zunbdb6 , ZUNBDB6 )
#define LAPACK_EXPORT_zunbdb F77_FUNC( zunbdb , ZUNBDB )
#define LAPACK_EXPORT_zuncsd2by1 F77_FUNC( zuncsd2by1 , ZUNCSD2BY1 )
#define LAPACK_EXPORT_zuncsd F77_FUNC( zuncsd , ZUNCSD )
#define LAPACK_EXPORT_zung2l F77_FUNC( zung2l , ZUNG2L )
#define LAPACK_EXPORT_zung2r F77_FUNC( zung2r , ZUNG2R )
#define LAPACK_EXPORT_zungbr F77_FUNC( zungbr , ZUNGBR )
#define LAPACK_EXPORT_zunghr F77_FUNC( zunghr , ZUNGHR )
#define LAPACK_EXPORT_zungl2 F77_FUNC( zungl2 , ZUNGL2 )
#define LAPACK_EXPORT_zunglq F77_FUNC( zunglq , ZUNGLQ )
#define LAPACK_EXPORT_zungql F77_FUNC( zungql , ZUNGQL )
#define LAPACK_EXPORT_zungqr F77_FUNC( zungqr , ZUNGQR )
#define LAPACK_EXPORT_zungr2 F77_FUNC( zungr2 , ZUNGR2 )
#define LAPACK_EXPORT_zungrq F77_FUNC( zungrq , ZUNGRQ )
#define LAPACK_EXPORT_zungtr F77_FUNC( zungtr , ZUNGTR )
#define LAPACK_EXPORT_zungtsqr F77_FUNC( zungtsqr , ZUNGTSQR )
#define LAPACK_EXPORT_zunhr_col F77_FUNC( zunhr_col , ZUNHR_COL )
#define LAPACK_EXPORT_zunm22 F77_FUNC( zunm22 , ZUNM22 )
#define LAPACK_EXPORT_zunm2l F77_FUNC( zunm2l , ZUNM2L )
#define LAPACK_EXPORT_zunm2r F77_FUNC( zunm2r , ZUNM2R )
#define LAPACK_EXPORT_zunmbr F77_FUNC( zunmbr , ZUNMBR )
#define LAPACK_EXPORT_zunmhr F77_FUNC( zunmhr , ZUNMHR )
#define LAPACK_EXPORT_zunml2 F77_FUNC( zunml2 , ZUNML2 )
#define LAPACK_EXPORT_zunmlq F77_FUNC( zunmlq , ZUNMLQ )
#define LAPACK_EXPORT_zunmql F77_FUNC( zunmql , ZUNMQL )
#define LAPACK_EXPORT_zunmqr F77_FUNC( zunmqr , ZUNMQR )
#define LAPACK_EXPORT_zunmr2 F77_FUNC( zunmr2 , ZUNMR2 )
#define LAPACK_EXPORT_zunmr3 F77_FUNC( zunmr3 , ZUNMR3 )
#define LAPACK_EXPORT_zunmrq F77_FUNC( zunmrq , ZUNMRQ )
#define LAPACK_EXPORT_zunmrz F77_FUNC( zunmrz , ZUNMRZ )
#define LAPACK_EXPORT_zunmtr F77_FUNC( zunmtr , ZUNMTR )
#define LAPACK_EXPORT_zupgtr F77_FUNC( zupgtr , ZUPGTR )
#define LAPACK_EXPORT_zupmtr F77_FUNC( zupmtr , ZUPMTR )
#define LAPACK_EXPORT_sgetrfnp F77_FUNC( sgetrfnp , SGETRFNP )
#define LAPACK_EXPORT_dgetrfnp F77_FUNC( dgetrfnp , DGETRFNP )
#define LAPACK_EXPORT_cgetrfnp F77_FUNC( cgetrfnp , CGETRFNP )
#define LAPACK_EXPORT_zgetrfnp F77_FUNC( zgetrfnp , ZGETRFNP )
#define LAPACK_EXPORT_sspffrt2 F77_FUNC( sspffrt2 , SSPFFRT2 )
#define LAPACK_EXPORT_dspffrt2 F77_FUNC( dspffrt2 , DSPFFRT2 )
#define LAPACK_EXPORT_cspffrt2 F77_FUNC( cspffrt2 , CSPFFRT2 )
#define LAPACK_EXPORT_zspffrt2 F77_FUNC( zspffrt2 , ZSPFFRT2 )
#define LAPACK_EXPORT_sspffrtx F77_FUNC( sspffrtx , SSPFFRTX )
#define LAPACK_EXPORT_dspffrtx F77_FUNC( dspffrtx , DSPFFRTX )
#define LAPACK_EXPORT_cspffrtx F77_FUNC( cspffrtx , CSPFFRTX )
#define LAPACK_EXPORT_zspffrtx F77_FUNC( zspffrtx , ZSPFFRTX )
#define LAPACK_EXPORT_cgetrfnpi F77_FUNC( cgetrfnpi , CGETRFNPI )
#define LAPACK_EXPORT_dgetrfnpi F77_FUNC( dgetrfnpi , DGETRFNPI )
#define LAPACK_EXPORT_sgetrfnpi F77_FUNC( sgetrfnpi , SGETRFNPI )
#define LAPACK_EXPORT_zgetrfnpi F77_FUNC( zgetrfnpi , ZGETRFNPI )

#define LAPACK_EXPORT_sopgtr F77_FUNC( sopgtr, SOPGTR )
#define LAPACK_EXPORT_dopgtr F77_FUNC( dopgtr, DOPGTR )
#define LAPACK_EXPORT_sorcsd F77_FUNC( sorcsd, SORCSD )
#define LAPACK_EXPORT_dorcsd F77_FUNC( dorcsd, DORCSD )
#define LAPACK_EXPORT_sorcsd2by1 F77_FUNC( sorcsd2by1, SORCSD2BY1 )
#define LAPACK_EXPORT_dorcsd2by1 F77_FUNC( dorcsd2by1, DORCSD2BY1 )
#define LAPACK_EXPORT_sorghr F77_FUNC( sorghr, SORGHR )
#define LAPACK_EXPORT_dorghr F77_FUNC( dorghr, DORGHR )
#define LAPACK_EXPORT_sormhr F77_FUNC( sormhr, SORMHR )
#define LAPACK_EXPORT_dormhr F77_FUNC( dormhr, DORMHR )

#define LAPACK_EXPORT_sgedmd F77_FUNC( sgedmd, SGEDMD )
#define LAPACK_EXPORT_dgedmd F77_FUNC( dgedmd, DGEDMD )
#define LAPACK_EXPORT_cgedmd F77_FUNC( cgedmd, CGEDMD )
#define LAPACK_EXPORT_zgedmd F77_FUNC( zgedmd, ZGEDMD )

#define LAPACK_EXPORT_sgedmdq F77_FUNC( sgedmdq, SGEDMDQ )
#define LAPACK_EXPORT_dgedmdq F77_FUNC( dgedmdq, DGEDMDQ )
#define LAPACK_EXPORT_cgedmdq F77_FUNC( cgedmdq, CGEDMDQ )
#define LAPACK_EXPORT_zgedmdq F77_FUNC( zgedmdq, ZGEDMDQ )




// Function Prototypes declaration
void LAPACK_EXPORT_cgelst(char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t * nrhs, scomplex *a, aocl_int_t *lda, scomplex *b, aocl_int_t *ldb, scomplex * work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_clatrs3(char *uplo, char *trans, char *diag, char * normin, aocl_int_t *n, aocl_int_t *nrhs, scomplex *a, aocl_int_t *lda, scomplex * x, aocl_int_t *ldx, real *scale, real *cnorm, real *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_ctrsyl3(char *trana, char *tranb, aocl_int_t *isgn, aocl_int_t *m, aocl_int_t *n, scomplex *a, aocl_int_t *lda, scomplex *b, aocl_int_t *ldb, scomplex *c__, aocl_int_t *ldc, real *scale, real *swork, aocl_int_t * ldswork, aocl_int_t *info);
void LAPACK_EXPORT_dgelst(char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t * nrhs, doublereal *a, aocl_int_t *lda, doublereal *b, aocl_int_t *ldb, doublereal *work, aocl_int_t *lwork, aocl_int_t *info);
doublereal LAPACK_EXPORT_dlarmm(doublereal *anorm, doublereal *bnorm, doublereal *cnorm);
void LAPACK_EXPORT_dlatrs3(char *uplo, char *trans, char *diag, char * normin, aocl_int_t *n, aocl_int_t *nrhs, doublereal *a, aocl_int_t *lda, doublereal *x, aocl_int_t *ldx, doublereal *scale, doublereal *cnorm, doublereal *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_dtrsyl3(char *trana, char *tranb, aocl_int_t *isgn, aocl_int_t *m, aocl_int_t *n, doublereal *a, aocl_int_t *lda, doublereal *b, aocl_int_t *ldb, doublereal *c__, aocl_int_t *ldc, doublereal *scale, aocl_int_t *iwork, aocl_int_t *liwork, doublereal *swork, aocl_int_t *ldswork, aocl_int_t *info);
void LAPACK_EXPORT_sgelst(char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t * nrhs, real *a, aocl_int_t *lda, real *b, aocl_int_t *ldb, real *work, aocl_int_t *lwork, aocl_int_t *info);
real LAPACK_EXPORT_slarmm(real *anorm, real *bnorm, real *cnorm);
void LAPACK_EXPORT_slatrs3(char *uplo, char *trans, char *diag, char * normin, aocl_int_t *n, aocl_int_t *nrhs, real *a, aocl_int_t *lda, real *x, aocl_int_t *ldx, real *scale, real *cnorm, real *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_strsyl3(char *trana, char *tranb, aocl_int_t *isgn, aocl_int_t *m, aocl_int_t *n, real *a, aocl_int_t *lda, real *b, aocl_int_t *ldb, real *c__, aocl_int_t *ldc, real *scale, aocl_int_t *iwork, aocl_int_t *liwork, real *swork, aocl_int_t *ldswork, aocl_int_t *info);
void LAPACK_EXPORT_zgelst(char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t * nrhs, dcomplex *a, aocl_int_t *lda, dcomplex *b, aocl_int_t *ldb, dcomplex *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_zlatrs3(char *uplo, char *trans, char *diag, char * normin, aocl_int_t *n, aocl_int_t *nrhs, dcomplex *a, aocl_int_t *lda, dcomplex *x, aocl_int_t *ldx, doublereal *scale, doublereal *cnorm, doublereal *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_ztrsyl3(char *trana, char *tranb, aocl_int_t *isgn, aocl_int_t *m, aocl_int_t *n, dcomplex *a, aocl_int_t *lda, dcomplex *b, aocl_int_t *ldb, dcomplex *c__, aocl_int_t *ldc, doublereal *scale, doublereal *swork, aocl_int_t *ldswork, aocl_int_t *info);
double LAPACK_EXPORT_dlamch(char* cmach);
double LAPACK_EXPORT_dlamc3(double* a, double* b);
doublereal LAPACK_EXPORT_dladiv2(doublereal* a, doublereal* b, doublereal* c__, doublereal* d__, doublereal* r__, doublereal* t);
void LAPACK_EXPORT_dladiv1(doublereal* a, doublereal* b, doublereal* c__, doublereal* d__, doublereal* p, doublereal* q);
real LAPACK_EXPORT_sladiv2(float *a, float* b, float* c, float* d, float* r, float* t);
void LAPACK_EXPORT_sladiv1(float* a, float* b, float* c, float* d, float* p, float* q);
aocl_int_t LAPACK_EXPORT_iparmq(aocl_int_t* ispec, char* name, char* opts, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, aocl_int_t* lwork);
aocl_int_t LAPACK_EXPORT_ilaenv(aocl_int_t* ispec, char* name, char* opts, aocl_int_t* n1, aocl_int_t* n2, aocl_int_t* n3, aocl_int_t* n4);
aocl_int_t LAPACK_EXPORT_ieeeck(aocl_int_t* ispec, float* zero, float* one);
logical LAPACK_EXPORT_lsamen(const aocl_int_t* n, const char* ca, const char* cb);
float LAPACK_EXPORT_slamc3(float* a, float* b);
float LAPACK_EXPORT_slamch(char* cmach);
void LAPACK_EXPORT_cgetsqrhrt(aocl_int_t *m, aocl_int_t *n, aocl_int_t *mb1, aocl_int_t *nb1, aocl_int_t *nb2, scomplex *a, aocl_int_t *lda, scomplex *t, aocl_int_t *ldt, scomplex *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_claqz0(char *wants, char *wantq, char *wantz, aocl_int_t * n, aocl_int_t *ilo, aocl_int_t *ihi, scomplex *a, aocl_int_t *lda, scomplex *b, aocl_int_t *ldb, scomplex *alpha, scomplex *beta, scomplex *q, aocl_int_t *ldq, scomplex *z__, aocl_int_t *ldz, scomplex *work, aocl_int_t *lwork, real * rwork, aocl_int_t *rec, aocl_int_t *info);
void LAPACK_EXPORT_claqz1(logical *ilq, logical *ilz, aocl_int_t *k, aocl_int_t * istartm, aocl_int_t *istopm, aocl_int_t *ihi, scomplex *a, aocl_int_t *lda, scomplex *b, aocl_int_t *ldb, aocl_int_t *nq, aocl_int_t *qstart, scomplex *q, aocl_int_t *ldq, aocl_int_t *nz, aocl_int_t *zstart, scomplex *z__, aocl_int_t * ldz);
void LAPACK_EXPORT_claqz2(logical *ilschur, logical *ilq, logical *ilz, aocl_int_t *n, aocl_int_t *ilo, aocl_int_t *ihi, aocl_int_t *nw, scomplex *a, aocl_int_t *lda, scomplex *b, aocl_int_t *ldb, scomplex *q, aocl_int_t *ldq, scomplex *z__, aocl_int_t *ldz, aocl_int_t *ns, aocl_int_t *nd, scomplex *alpha, scomplex *beta, scomplex *qc, aocl_int_t *ldqc, scomplex *zc, aocl_int_t *ldzc, scomplex *work, aocl_int_t *lwork, real *rwork, aocl_int_t *rec, aocl_int_t * info);
void LAPACK_EXPORT_claqz3(logical *ilschur, logical *ilq, logical *ilz, aocl_int_t *n, aocl_int_t *ilo, aocl_int_t *ihi, aocl_int_t *nshifts, aocl_int_t * nblock_desired__, scomplex *alpha, scomplex *beta, scomplex *a, aocl_int_t * lda, scomplex *b, aocl_int_t *ldb, scomplex *q, aocl_int_t *ldq, scomplex *z__, aocl_int_t *ldz, scomplex *qc, aocl_int_t *ldqc, scomplex *zc, aocl_int_t *ldzc, scomplex *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_clarfb_gett(char *ident, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, scomplex *t, aocl_int_t *ldt, scomplex *a, aocl_int_t *lda, scomplex *b, aocl_int_t *ldb, scomplex *work, aocl_int_t *ldwork);
void LAPACK_EXPORT_cungtsqr_row(aocl_int_t *m, aocl_int_t *n, aocl_int_t *mb, aocl_int_t *nb, scomplex *a, aocl_int_t *lda, scomplex *t, aocl_int_t *ldt, scomplex *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_dgetsqrhrt(aocl_int_t *m, aocl_int_t *n, aocl_int_t *mb1, aocl_int_t *nb1, aocl_int_t *nb2, doublereal *a, aocl_int_t *lda, doublereal * t, aocl_int_t *ldt, doublereal *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_dlaqz0(char *wants, char *wantq, char *wantz, aocl_int_t * n, aocl_int_t *ilo, aocl_int_t *ihi, doublereal *a, aocl_int_t *lda, doublereal *b, aocl_int_t *ldb, doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal *q, aocl_int_t *ldq, doublereal *z__, aocl_int_t *ldz, doublereal *work, aocl_int_t *lwork, aocl_int_t *rec, aocl_int_t *info);
void LAPACK_EXPORT_dlaqz1(doublereal *a, aocl_int_t *lda, doublereal *b, aocl_int_t *ldb, doublereal *sr1, doublereal *sr2, doublereal *si, doublereal *beta1, doublereal *beta2, doublereal *v);
void LAPACK_EXPORT_dlaqz2(logical *ilq, logical *ilz, aocl_int_t *k, aocl_int_t * istartm, aocl_int_t *istopm, aocl_int_t *ihi, doublereal *a, aocl_int_t *lda, doublereal *b, aocl_int_t *ldb, aocl_int_t *nq, aocl_int_t *qstart, doublereal *q, aocl_int_t *ldq, aocl_int_t *nz, aocl_int_t *zstart, doublereal *z__, aocl_int_t *ldz);
void LAPACK_EXPORT_dlaqz3(logical *ilschur, logical *ilq, logical *ilz, aocl_int_t *n, aocl_int_t *ilo, aocl_int_t *ihi, aocl_int_t *nw, doublereal *a, aocl_int_t *lda, doublereal *b, aocl_int_t *ldb, doublereal *q, aocl_int_t * ldq, doublereal *z__, aocl_int_t *ldz, aocl_int_t *ns, aocl_int_t *nd, doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal * qc, aocl_int_t *ldqc, doublereal *zc, aocl_int_t *ldzc, doublereal *work, aocl_int_t *lwork, aocl_int_t *rec, aocl_int_t *info);
void LAPACK_EXPORT_dlaqz4(logical *ilschur, logical *ilq, logical *ilz, aocl_int_t *n, aocl_int_t *ilo, aocl_int_t *ihi, aocl_int_t *nshifts, aocl_int_t * nblock_desired__, doublereal *sr, doublereal *si, doublereal *ss, doublereal *a, aocl_int_t *lda, doublereal *b, aocl_int_t *ldb, doublereal * q, aocl_int_t *ldq, doublereal *z__, aocl_int_t *ldz, doublereal *qc, aocl_int_t *ldqc, doublereal *zc, aocl_int_t *ldzc, doublereal *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_dlarfb_gett(char *ident, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, doublereal *t, aocl_int_t *ldt, doublereal *a, aocl_int_t *lda, doublereal *b, aocl_int_t *ldb, doublereal *work, aocl_int_t *ldwork);
void LAPACK_EXPORT_dorgtsqr_row(aocl_int_t *m, aocl_int_t *n, aocl_int_t *mb, aocl_int_t *nb, doublereal *a, aocl_int_t *lda, doublereal *t, aocl_int_t *ldt, doublereal *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_sgetsqrhrt(aocl_int_t *m, aocl_int_t *n, aocl_int_t *mb1, aocl_int_t *nb1, aocl_int_t *nb2, real *a, aocl_int_t *lda, real *t, aocl_int_t * ldt, real *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_slaqz0(char *wants, char *wantq, char *wantz, aocl_int_t * n, aocl_int_t *ilo, aocl_int_t *ihi, real *a, aocl_int_t *lda, real *b, aocl_int_t *ldb, real *alphar, real *alphai, real *beta, real *q, aocl_int_t *ldq, real *z__, aocl_int_t *ldz, real *work, aocl_int_t *lwork, aocl_int_t *rec, aocl_int_t *info);
void LAPACK_EXPORT_slaqz1(real *a, aocl_int_t *lda, real *b, aocl_int_t *ldb, real *sr1, real *sr2, real *si, real *beta1, real *beta2, real *v);
void LAPACK_EXPORT_slaqz2(logical *ilq, logical *ilz, aocl_int_t *k, aocl_int_t * istartm, aocl_int_t *istopm, aocl_int_t *ihi, real *a, aocl_int_t *lda, real * b, aocl_int_t *ldb, aocl_int_t *nq, aocl_int_t *qstart, real *q, aocl_int_t *ldq, aocl_int_t *nz, aocl_int_t *zstart, real *z__, aocl_int_t *ldz);
void LAPACK_EXPORT_slaqz3(logical *ilschur, logical *ilq, logical *ilz, aocl_int_t *n, aocl_int_t *ilo, aocl_int_t *ihi, aocl_int_t *nw, real *a, aocl_int_t *lda, real *b, aocl_int_t *ldb, real *q, aocl_int_t *ldq, real *z__, aocl_int_t *ldz, aocl_int_t *ns, aocl_int_t *nd, real *alphar, real *alphai, real *beta, real *qc, aocl_int_t *ldqc, real *zc, aocl_int_t *ldzc, real * work, aocl_int_t *lwork, aocl_int_t *rec, aocl_int_t *info);
void LAPACK_EXPORT_slaqz4(logical *ilschur, logical *ilq, logical *ilz, aocl_int_t *n, aocl_int_t *ilo, aocl_int_t *ihi, aocl_int_t *nshifts, aocl_int_t * nblock_desired__, real *sr, real *si, real *ss, real *a, aocl_int_t *lda, real *b, aocl_int_t *ldb, real *q, aocl_int_t *ldq, real *z__, aocl_int_t * ldz, real *qc, aocl_int_t *ldqc, real *zc, aocl_int_t *ldzc, real *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_slarfb_gett(char *ident, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *t, aocl_int_t *ldt, real *a, aocl_int_t *lda, real *b, aocl_int_t *ldb, real *work, aocl_int_t *ldwork);
void LAPACK_EXPORT_sorgtsqr_row(aocl_int_t *m, aocl_int_t *n, aocl_int_t *mb, aocl_int_t *nb, real *a, aocl_int_t *lda, real *t, aocl_int_t *ldt, real *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_zgetsqrhrt(aocl_int_t *m, aocl_int_t *n, aocl_int_t *mb1, aocl_int_t *nb1, aocl_int_t *nb2, dcomplex *a, aocl_int_t *lda, dcomplex *t, aocl_int_t *ldt, dcomplex *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_zlaqz0(char *wants, char *wantq, char *wantz, aocl_int_t * n, aocl_int_t *ilo, aocl_int_t *ihi, dcomplex *a, aocl_int_t *lda, dcomplex *b, aocl_int_t *ldb, dcomplex *alpha, dcomplex * beta, dcomplex *q, aocl_int_t *ldq, dcomplex *z__, aocl_int_t * ldz, dcomplex *work, aocl_int_t *lwork, doublereal *rwork, aocl_int_t * rec, aocl_int_t *info);
void LAPACK_EXPORT_zlaqz1(logical *ilq, logical *ilz, aocl_int_t *k, aocl_int_t * istartm, aocl_int_t *istopm, aocl_int_t *ihi, dcomplex *a, aocl_int_t * lda, dcomplex *b, aocl_int_t *ldb, aocl_int_t *nq, aocl_int_t *qstart, dcomplex *q, aocl_int_t *ldq, aocl_int_t *nz, aocl_int_t *zstart, dcomplex *z__, aocl_int_t *ldz);
void LAPACK_EXPORT_zlaqz2(logical *ilschur, logical *ilq, logical *ilz, aocl_int_t *n, aocl_int_t *ilo, aocl_int_t *ihi, aocl_int_t *nw, dcomplex *a, aocl_int_t *lda, dcomplex *b, aocl_int_t *ldb, dcomplex *q, aocl_int_t *ldq, dcomplex *z__, aocl_int_t *ldz, aocl_int_t *ns, aocl_int_t * nd, dcomplex *alpha, dcomplex *beta, dcomplex *qc, aocl_int_t *ldqc, dcomplex *zc, aocl_int_t *ldzc, dcomplex *work, aocl_int_t *lwork, doublereal *rwork, aocl_int_t *rec, aocl_int_t *info);
void LAPACK_EXPORT_zlaqz3(logical *ilschur, logical *ilq, logical *ilz, aocl_int_t *n, aocl_int_t *ilo, aocl_int_t *ihi, aocl_int_t *nshifts, aocl_int_t * nblock_desired__, dcomplex *alpha, dcomplex *beta, dcomplex *a, aocl_int_t *lda, dcomplex *b, aocl_int_t *ldb, dcomplex *q, aocl_int_t *ldq, dcomplex *z__, aocl_int_t *ldz, dcomplex *qc, aocl_int_t *ldqc, dcomplex *zc, aocl_int_t *ldzc, dcomplex *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_zlarfb_gett(char *ident, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, dcomplex *t, aocl_int_t *ldt, dcomplex *a, aocl_int_t *lda, dcomplex *b, aocl_int_t *ldb, dcomplex *work, aocl_int_t * ldwork);
void LAPACK_EXPORT_zungtsqr_row(aocl_int_t *m, aocl_int_t *n, aocl_int_t *mb, aocl_int_t *nb, dcomplex *a, aocl_int_t *lda, dcomplex *t, aocl_int_t *ldt, dcomplex *work, aocl_int_t *lwork, aocl_int_t *info);
double LAPACK_EXPORT_dla_gbrcond(char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, double* ab, aocl_int_t* ldab, double* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, aocl_int_t* cmode, double* c, aocl_int_t* info, double* work, aocl_int_t* iwork);
double LAPACK_EXPORT_dla_gbrpvgrw(aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* ncols, double* ab, aocl_int_t* ldab, double* afb, aocl_int_t* ldafb);
double LAPACK_EXPORT_dla_gercond(char* trans, aocl_int_t* n, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, aocl_int_t* ipiv, aocl_int_t* cmode, double* c, aocl_int_t* info, double* work, aocl_int_t* iwork);
double LAPACK_EXPORT_dla_gerpvgrw(aocl_int_t* n, aocl_int_t* ncols, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf);
double LAPACK_EXPORT_dlangb(char* norm, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, double* ab, aocl_int_t* ldab, double* work);
double LAPACK_EXPORT_dlange(char* norm, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* work);
double LAPACK_EXPORT_dlangt(char* norm, aocl_int_t* n, double* dl, double* d, double* du);
double LAPACK_EXPORT_dlanhs(char* norm, aocl_int_t* n, double* a, aocl_int_t* lda, double* work);
double LAPACK_EXPORT_dlansb(char* norm, char* uplo, aocl_int_t* n, aocl_int_t* k, double* ab, aocl_int_t* ldab, double* work);
double LAPACK_EXPORT_dlansf(char* norm, char* transr, char* uplo, aocl_int_t* n, double* a, double* work);
double LAPACK_EXPORT_dlansp(char* norm, char* uplo, aocl_int_t* n, double* ap, double* work);
double LAPACK_EXPORT_dlanst(char* norm, aocl_int_t* n, double* d, double* e);
double LAPACK_EXPORT_dlansy(char* norm, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* work);
double LAPACK_EXPORT_dlantb(char* norm, char* uplo, char* diag, aocl_int_t* n, aocl_int_t* k, double* ab, aocl_int_t* ldab, double* work);
double LAPACK_EXPORT_dlantp(char* norm, char* uplo, char* diag, aocl_int_t* n, double* ap, double* work);
double LAPACK_EXPORT_dlantr(char* norm, char* uplo, char* diag, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* work);
double LAPACK_EXPORT_dla_porcond(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, aocl_int_t* cmode, double* c, aocl_int_t* info, double* work, aocl_int_t* iwork);
double LAPACK_EXPORT_dla_porpvgrw(char* uplo, aocl_int_t* ncols, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, double* work);
double LAPACK_EXPORT_dlapy2(double* x, double* y);
double LAPACK_EXPORT_dlapy3(double* x, double* y, double* z__);
double LAPACK_EXPORT_dla_syrcond(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, aocl_int_t* ipiv, aocl_int_t* cmode, double* c, aocl_int_t* info, double* work, aocl_int_t* iwork);
double LAPACK_EXPORT_dla_syrpvgrw(char* uplo, aocl_int_t* n, aocl_int_t* info, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, aocl_int_t* ipiv, double* work);
double LAPACK_EXPORT_dzsum1(aocl_int_t* n, dcomplex* cx, aocl_int_t* incx);
double LAPACK_EXPORT_zla_gbrcond_c(char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, dcomplex* ab, aocl_int_t* ldab, dcomplex* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, double* c, logical* capply, aocl_int_t* info, dcomplex* work, double* rwork);
double LAPACK_EXPORT_zla_gbrcond_x(char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, dcomplex* ab, aocl_int_t* ldab, dcomplex* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, dcomplex* x, aocl_int_t* info, dcomplex* work, double* rwork);
double LAPACK_EXPORT_zla_gbrpvgrw(aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* ncols, dcomplex* ab, aocl_int_t* ldab, dcomplex* afb, aocl_int_t* ldafb);
double LAPACK_EXPORT_zla_gercond_c(char* trans, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, double* c, logical* capply, aocl_int_t* info, dcomplex* work, double* rwork);
double LAPACK_EXPORT_zla_gercond_x(char* trans, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, dcomplex* x, aocl_int_t* info, dcomplex* work, double* rwork);
double LAPACK_EXPORT_zla_gerpvgrw(aocl_int_t* n, aocl_int_t* ncols, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf);
double LAPACK_EXPORT_zla_hercond_c(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, double* c, logical* capply, aocl_int_t* info, dcomplex* work, double* rwork);
double LAPACK_EXPORT_zla_hercond_x(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, dcomplex* x, aocl_int_t* info, dcomplex* work, double* rwork);
double LAPACK_EXPORT_zla_herpvgrw(char* uplo, aocl_int_t* n, aocl_int_t* info, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, double* work);
double LAPACK_EXPORT_zlangb(char* norm, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, dcomplex* ab, aocl_int_t* ldab, double* work);
double LAPACK_EXPORT_zlange(char* norm, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* work);
double LAPACK_EXPORT_zlangt(char* norm, aocl_int_t* n, dcomplex* dl, dcomplex* d, dcomplex* du);
double LAPACK_EXPORT_zlanhb(char* norm, char* uplo, aocl_int_t* n, aocl_int_t* k, dcomplex* ab, aocl_int_t* ldab, double* work);
double LAPACK_EXPORT_zlanhe(char* norm, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* work);
double LAPACK_EXPORT_zlanhf(char* norm, char* transr, char* uplo, aocl_int_t* n, dcomplex* a, double* work);
double LAPACK_EXPORT_zlanhp(char* norm, char* uplo, aocl_int_t* n, dcomplex* ap, double* work);
double LAPACK_EXPORT_zlanhs(char* norm, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* work);
double LAPACK_EXPORT_zlanht(char* norm, aocl_int_t* n, double* d, dcomplex* e);
double LAPACK_EXPORT_zlansb(char* norm, char* uplo, aocl_int_t* n, aocl_int_t* k, dcomplex* ab, aocl_int_t* ldab, double* work);
double LAPACK_EXPORT_zlansp(char* norm, char* uplo, aocl_int_t* n, dcomplex* ap, double* work);
double LAPACK_EXPORT_zlansy(char* norm, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* work);
double LAPACK_EXPORT_zlantb(char* norm, char* uplo, char* diag, aocl_int_t* n, aocl_int_t* k, dcomplex* ab, aocl_int_t* ldab, double* work);
double LAPACK_EXPORT_zlantp(char* norm, char* uplo, char* diag, aocl_int_t* n, dcomplex* ap, double* work);
double LAPACK_EXPORT_zlantr(char* norm, char* uplo, char* diag, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* work);
double LAPACK_EXPORT_zla_porcond_c(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, double* c, logical* capply, aocl_int_t* info, dcomplex* work, double* rwork);
double LAPACK_EXPORT_zla_porcond_x(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, dcomplex* x, aocl_int_t* info, dcomplex* work, double* rwork);
double LAPACK_EXPORT_zla_porpvgrw(char* uplo, aocl_int_t* ncols, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, double* work);
double LAPACK_EXPORT_zla_syrcond_c(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, double* c, logical* capply, aocl_int_t* info, dcomplex* work, double* rwork);
double LAPACK_EXPORT_zla_syrcond_x(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, dcomplex* x, aocl_int_t* info, dcomplex* work, double* rwork);
double LAPACK_EXPORT_zla_syrpvgrw(char* uplo, aocl_int_t* n, aocl_int_t* info, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, double* work);
float LAPACK_EXPORT_cla_gbrcond_c(char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, scomplex* ab, aocl_int_t* ldab, scomplex* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, float* c, logical* capply, aocl_int_t* info, scomplex* work, float* rwork);
float LAPACK_EXPORT_cla_gbrcond_x(char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, scomplex* ab, aocl_int_t* ldab, scomplex* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, scomplex* x, aocl_int_t* info, scomplex* work, float* rwork);
float LAPACK_EXPORT_cla_gbrpvgrw(aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* ncols, scomplex* ab, aocl_int_t* ldab, scomplex* afb, aocl_int_t* ldafb);
float LAPACK_EXPORT_cla_gercond_c(char* trans, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, float* c, logical* capply, aocl_int_t* info, scomplex* work, float* rwork);
float LAPACK_EXPORT_cla_gercond_x(char* trans, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, scomplex* x, aocl_int_t* info, scomplex* work, float* rwork);
float LAPACK_EXPORT_cla_gerpvgrw(aocl_int_t* n, aocl_int_t* ncols, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf);
float LAPACK_EXPORT_cla_hercond_c(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, float* c, logical* capply, aocl_int_t* info, scomplex* work, float* rwork);
float LAPACK_EXPORT_cla_hercond_x(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, scomplex* x, aocl_int_t* info, scomplex* work, float* rwork);
float LAPACK_EXPORT_cla_herpvgrw(char* uplo, aocl_int_t* n, aocl_int_t* info, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, float* work);
float LAPACK_EXPORT_clangb(char* norm, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, scomplex* ab, aocl_int_t* ldab, float* work);
float LAPACK_EXPORT_clange(char* norm, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* work);
float LAPACK_EXPORT_clangt(char* norm, aocl_int_t* n, scomplex* dl, scomplex* d, scomplex* du);
float LAPACK_EXPORT_clanhb(char* norm, char* uplo, aocl_int_t* n, aocl_int_t* k, scomplex* ab, aocl_int_t* ldab, float* work);
float LAPACK_EXPORT_clanhe(char* norm, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* work);
float LAPACK_EXPORT_clanhf(char* norm, char* transr, char* uplo, aocl_int_t* n, scomplex* a, float* work);
float LAPACK_EXPORT_clanhp(char* norm, char* uplo, aocl_int_t* n, scomplex* ap, float* work);
float LAPACK_EXPORT_clanhs(char* norm, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* work);
float LAPACK_EXPORT_clanht(char* norm, aocl_int_t* n, float* d, scomplex* e);
float LAPACK_EXPORT_clansb(char* norm, char* uplo, aocl_int_t* n, aocl_int_t* k, scomplex* ab, aocl_int_t* ldab, float* work);
float LAPACK_EXPORT_clansp(char* norm, char* uplo, aocl_int_t* n, scomplex* ap, float* work);
float LAPACK_EXPORT_clansy(char* norm, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* work);
float LAPACK_EXPORT_clantb(char* norm, char* uplo, char* diag, aocl_int_t* n, aocl_int_t* k, scomplex* ab, aocl_int_t* ldab, float* work);
float LAPACK_EXPORT_clantp(char* norm, char* uplo, char* diag, aocl_int_t* n, scomplex* ap, float* work);
float LAPACK_EXPORT_clantr(char* norm, char* uplo, char* diag, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* work);
float LAPACK_EXPORT_cla_porcond_c(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, float* c, logical* capply, aocl_int_t* info, scomplex* work, float* rwork);
float LAPACK_EXPORT_cla_porcond_x(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, scomplex* x, aocl_int_t* info, scomplex* work, float* rwork);
float LAPACK_EXPORT_cla_porpvgrw(char* uplo, aocl_int_t* ncols, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, float* work);
float LAPACK_EXPORT_cla_syrcond_c(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, float* c, logical* capply, aocl_int_t* info, scomplex* work, float* rwork);
float LAPACK_EXPORT_cla_syrcond_x(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, scomplex* x, aocl_int_t* info, scomplex* work, float* rwork);
float LAPACK_EXPORT_cla_syrpvgrw(char* uplo, aocl_int_t* n, aocl_int_t* info, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, float* work);
float LAPACK_EXPORT_scsum1(aocl_int_t* n, scomplex* cx, aocl_int_t* incx);
float LAPACK_EXPORT_sla_gbrcond(char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, float* ab, aocl_int_t* ldab, float* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, aocl_int_t* cmode, float* c, aocl_int_t* info, float* work, aocl_int_t* iwork);
float LAPACK_EXPORT_sla_gbrpvgrw(aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* ncols, float* ab, aocl_int_t* ldab, float* afb, aocl_int_t* ldafb);
float LAPACK_EXPORT_sla_gercond(char* trans, aocl_int_t* n, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, aocl_int_t* ipiv, aocl_int_t* cmode, float* c, aocl_int_t* info, float* work, aocl_int_t* iwork);
float LAPACK_EXPORT_sla_gerpvgrw(aocl_int_t* n, aocl_int_t* ncols, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf);
float LAPACK_EXPORT_slangb(char* norm, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, float* ab, aocl_int_t* ldab, float* work);
float LAPACK_EXPORT_slange(char* norm, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* work);
float LAPACK_EXPORT_slangt(char* norm, aocl_int_t* n, float* dl, float* d, float* du);
float LAPACK_EXPORT_slanhs(char* norm, aocl_int_t* n, float* a, aocl_int_t* lda, float* work);
float LAPACK_EXPORT_slansb(char* norm, char* uplo, aocl_int_t* n, aocl_int_t* k, float* ab, aocl_int_t* ldab, float* work);
float LAPACK_EXPORT_slansf(char* norm, char* transr, char* uplo, aocl_int_t* n, float* a, float* work);
float LAPACK_EXPORT_slansp(char* norm, char* uplo, aocl_int_t* n, float* ap, float* work);
float LAPACK_EXPORT_slanst(char* norm, aocl_int_t* n, float* d, float* e);
float LAPACK_EXPORT_slansy(char* norm, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* work);
float LAPACK_EXPORT_slantb(char* norm, char* uplo, char* diag, aocl_int_t* n, aocl_int_t* k, float* ab, aocl_int_t* ldab, float* work);
float LAPACK_EXPORT_slantp(char* norm, char* uplo, char* diag, aocl_int_t* n, float* ap, float* work);
float LAPACK_EXPORT_slantr(char* norm, char* uplo, char* diag, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* work);
float LAPACK_EXPORT_sla_porcond(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, aocl_int_t* cmode, float* c, aocl_int_t* info, float* work, aocl_int_t* iwork);
float LAPACK_EXPORT_sla_porpvgrw(char* uplo, aocl_int_t* ncols, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, float* work);
float LAPACK_EXPORT_slapy2(float* x, float* y);
float LAPACK_EXPORT_slapy3(float* x, float* y, float* z__);
float LAPACK_EXPORT_sla_syrcond(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, aocl_int_t* ipiv, aocl_int_t* cmode, float* c, aocl_int_t* info, float* work, aocl_int_t* iwork);
float LAPACK_EXPORT_sla_syrpvgrw(char* uplo, aocl_int_t* n, aocl_int_t* info, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, aocl_int_t* ipiv, float* work);
void LAPACK_EXPORT_cbbcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, float* theta, float* phi, scomplex* u1, aocl_int_t* ldu1, scomplex* u2, aocl_int_t* ldu2, scomplex* v1t, aocl_int_t* ldv1t, scomplex* v2t, aocl_int_t* ldv2t, float* b11d, float* b11e, float* b12d, float* b12e, float* b21d, float* b21e, float* b22d, float* b22e, float* rwork, aocl_int_t* lrwork, aocl_int_t* info);
void LAPACK_EXPORT_cbdsqr(char* uplo, aocl_int_t* n, aocl_int_t* ncvt, aocl_int_t* nru, aocl_int_t* ncc, float* d, float* e, scomplex* vt, aocl_int_t* ldvt, scomplex* u, aocl_int_t* ldu, scomplex* c, aocl_int_t* ldc, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgbbrd(char* vect, aocl_int_t* m, aocl_int_t* n, aocl_int_t* ncc, aocl_int_t* kl, aocl_int_t* ku, scomplex* ab, aocl_int_t* ldab, float* d, float* e, scomplex* q, aocl_int_t* ldq, scomplex* pt, aocl_int_t* ldpt, scomplex* c, aocl_int_t* ldc, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgbcon(char* norm, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, scomplex* ab, aocl_int_t* ldab, aocl_int_t* ipiv, float* anorm, float* rcond, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgbequb(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, scomplex* ab, aocl_int_t* ldab, float* r, float* c, float* rowcnd, float* colcnd, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_cgbequ(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, scomplex* ab, aocl_int_t* ldab, float* r, float* c, float* rowcnd, float* colcnd, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_cgbrfs(char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, scomplex* ab, aocl_int_t* ldab, scomplex* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgbrfsx(char* trans, char* equed, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, scomplex* ab, aocl_int_t* ldab, scomplex* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, float* r, float* c, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgbsv(aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, scomplex* ab, aocl_int_t* ldab, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cgbsvx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, scomplex* ab, aocl_int_t* ldab, scomplex* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, char* equed, float* r, float* c, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgbsvxx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, scomplex* ab, aocl_int_t* ldab, scomplex* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, char* equed, float* r, float* c, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* rpvgrw, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgbtf2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, scomplex* ab, aocl_int_t* ldab, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_cgbtrf(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, scomplex* ab, aocl_int_t* ldab, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_cgbtrs(char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, scomplex* ab, aocl_int_t* ldab, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cgebak(char* job, char* side, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, float* scale, aocl_int_t* m, scomplex* v, aocl_int_t* ldv, aocl_int_t* info);
void LAPACK_EXPORT_cgebal(char* job, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ilo, aocl_int_t* ihi, float* scale, aocl_int_t* info);
void LAPACK_EXPORT_cgebd2(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* d, float* e, scomplex* tauq, scomplex* taup, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cgebrd(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* d, float* e, scomplex* tauq, scomplex* taup, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgecon(char* norm, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* anorm, float* rcond, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgeequb(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_cgeequ(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_cgees(char* jobvs, char* sort, L_fp select, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* sdim, scomplex* w, scomplex* vs, aocl_int_t* ldvs, scomplex* work, aocl_int_t* lwork, float* rwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_cgeesx(char* jobvs, char* sort, L_fp select, char* sense, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* sdim, scomplex* w, scomplex* vs, aocl_int_t* ldvs, float* rconde, float* rcondv, scomplex* work, aocl_int_t* lwork, float* rwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_cgeev(char* jobvl, char* jobvr, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* w, scomplex* vl, aocl_int_t* ldvl, scomplex* vr, aocl_int_t* ldvr, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgeevx(char* balanc, char* jobvl, char* jobvr, char* sense, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* w, scomplex* vl, aocl_int_t* ldvl, scomplex* vr, aocl_int_t* ldvr, aocl_int_t* ilo, aocl_int_t* ihi, float* scale, float* abnrm, float* rconde, float* rcondv, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgegs(char* jobvsl, char* jobvsr, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* alpha, scomplex* beta, scomplex* vsl, aocl_int_t* ldvsl, scomplex* vsr, aocl_int_t* ldvsr, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgegv(char* jobvl, char* jobvr, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* alpha, scomplex* beta, scomplex* vl, aocl_int_t* ldvl, scomplex* vr, aocl_int_t* ldvr, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgehd2(aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cgehrd(aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgejsv(char* joba, char* jobu, char* jobv, char* jobr, char* jobt, char* jobp, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* sva, scomplex* u, aocl_int_t* ldu, scomplex* v, aocl_int_t* ldv, scomplex* cwork, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_cgelq2(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cgelqf(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgelq(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* tsize, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgelqt3(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_cgelqt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* ldt, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cgels(char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgelsd(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, float* s, float* rcond, aocl_int_t* rank, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_cgelss(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, float* s, float* rcond, aocl_int_t* rank, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgelsx(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, aocl_int_t* jpvt, float* rcond, aocl_int_t* rank, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgelsy(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, aocl_int_t* jpvt, float* rcond, aocl_int_t* rank, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgemlq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* tsize, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgemlqt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* mb, scomplex* v, aocl_int_t* ldv, scomplex* t, aocl_int_t* ldt, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cgemqr(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* tsize, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgemqrt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* nb, scomplex* v, aocl_int_t* ldv, scomplex* t, aocl_int_t* ldt, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cgeql2(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cgeqlf(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgeqp3(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* jpvt, scomplex* tau, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgeqpf(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* jpvt, scomplex* tau, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgeqr2(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cgeqr2p(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cgeqrf(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgeqrfp(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgeqr(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* tsize, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgeqrt2(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_cgeqrt3(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_cgeqrt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nb, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* ldt, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cgerfs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgerfsx(char* trans, char* equed, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, float* r, float* c, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgerq2(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cgerqf(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgesc2(aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* rhs, aocl_int_t* ipiv, aocl_int_t* jpiv, float* scale);
void LAPACK_EXPORT_cgesdd(char* jobz, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* s, scomplex* u, aocl_int_t* ldu, scomplex* vt, aocl_int_t* ldvt, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_cgesvd(char* jobu, char* jobv, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* s, scomplex* u, aocl_int_t* ldu, scomplex* vt, aocl_int_t* ldvt, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgesvdq(char* joba, char* jobp, char* jobr, char* jobu, char* jobv, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* s, scomplex* u, aocl_int_t* ldu, scomplex* v, aocl_int_t* ldv, aocl_int_t* numrank, aocl_int_t* iwork, aocl_int_t* liwork, scomplex* cwork, aocl_int_t* lcwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* info);
void LAPACK_EXPORT_cgesvdx(char* jobu, char* jobvt, char* range, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, aocl_int_t* ns, float* s, scomplex* u, aocl_int_t* ldu, scomplex* vt, aocl_int_t* ldvt, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_cgesv(aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cgesvj(char* joba, char* jobu, char* jobv, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* sva, aocl_int_t* mv, scomplex* v, aocl_int_t* ldv, scomplex* cwork, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* info);
void LAPACK_EXPORT_cgesvx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, char* equed, float* r, float* c, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgesvxx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, char* equed, float* r, float* c, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* rpvgrw, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgetc2(aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* jpiv, aocl_int_t* info);
void LAPACK_EXPORT_cgetf2(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_cgetrf2(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_cgetrf(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_cgetri(aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgetrs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cgetsls(char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cggbak(char* job, char* side, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, float* lscale, float* rscale, aocl_int_t* m, scomplex* v, aocl_int_t* ldv, aocl_int_t* info);
void LAPACK_EXPORT_cggbal(char* job, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, aocl_int_t* ilo, aocl_int_t* ihi, float* lscale, float* rscale, float* work, aocl_int_t* info);
void LAPACK_EXPORT_cgges3(char* jobvsl, char* jobvsr, char* sort, L_fp selctg, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, aocl_int_t* sdim, scomplex* alpha, scomplex* beta, scomplex* vsl, aocl_int_t* ldvsl, scomplex* vsr, aocl_int_t* ldvsr, scomplex* work, aocl_int_t* lwork, float* rwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_cgges(char* jobvsl, char* jobvsr, char* sort, L_fp selctg, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, aocl_int_t* sdim, scomplex* alpha, scomplex* beta, scomplex* vsl, aocl_int_t* ldvsl, scomplex* vsr, aocl_int_t* ldvsr, scomplex* work, aocl_int_t* lwork, float* rwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_cggesx(char* jobvsl, char* jobvsr, char* sort, L_fp selctg, char* sense, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, aocl_int_t* sdim, scomplex* alpha, scomplex* beta, scomplex* vsl, aocl_int_t* ldvsl, scomplex* vsr, aocl_int_t* ldvsr, float* rconde, float* rcondv, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* iwork, aocl_int_t* liwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_cggev3(char* jobvl, char* jobvr, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* alpha, scomplex* beta, scomplex* vl, aocl_int_t* ldvl, scomplex* vr, aocl_int_t* ldvr, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cggev(char* jobvl, char* jobvr, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* alpha, scomplex* beta, scomplex* vl, aocl_int_t* ldvl, scomplex* vr, aocl_int_t* ldvr, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cggevx(char* balanc, char* jobvl, char* jobvr, char* sense, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* alpha, scomplex* beta, scomplex* vl, aocl_int_t* ldvl, scomplex* vr, aocl_int_t* ldvr, aocl_int_t* ilo, aocl_int_t* ihi, float* lscale, float* rscale, float* abnrm, float* bbnrm, float* rconde, float* rcondv, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* iwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_cggglm(aocl_int_t* n, aocl_int_t* m, aocl_int_t* p, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* d, scomplex* x, scomplex* y, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgghd3(char* compq, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* q, aocl_int_t* ldq, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgghrd(char* compq, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* q, aocl_int_t* ldq, scomplex* z, aocl_int_t* ldz, aocl_int_t* info);
void LAPACK_EXPORT_cgglse(aocl_int_t* m, aocl_int_t* n, aocl_int_t* p, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* c, scomplex* d, scomplex* x, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cggqrf(aocl_int_t* n, aocl_int_t* m, aocl_int_t* p, scomplex* a, aocl_int_t* lda, scomplex* taua, scomplex* b, aocl_int_t* ldb, scomplex* taub, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cggrqf(aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* taua, scomplex* b, aocl_int_t* ldb, scomplex* taub, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cggsvd3(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* n, aocl_int_t* p, aocl_int_t* k, aocl_int_t* l, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, float* alpha, float* beta, scomplex* u, aocl_int_t* ldu, scomplex* v, aocl_int_t* ldv, scomplex* q, aocl_int_t* ldq, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_cggsvd(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* n, aocl_int_t* p, aocl_int_t* k, aocl_int_t* l, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, float* alpha, float* beta, scomplex* u, aocl_int_t* ldu, scomplex* v, aocl_int_t* ldv, scomplex* q, aocl_int_t* ldq, scomplex* work, float* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_cggsvp3(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, float* tola, float* tolb, aocl_int_t* k, aocl_int_t* l, scomplex* u, aocl_int_t* ldu, scomplex* v, aocl_int_t* ldv, scomplex* q, aocl_int_t* ldq, aocl_int_t* iwork, float* rwork, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cggsvp(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, float* tola, float* tolb, aocl_int_t* k, aocl_int_t* l, scomplex* u, aocl_int_t* ldu, scomplex* v, aocl_int_t* ldv, scomplex* q, aocl_int_t* ldq, aocl_int_t* iwork, float* rwork, scomplex* tau, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cgsvj0(char* jobv, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* d, float* sva, aocl_int_t* mv, scomplex* v, aocl_int_t* ldv, float* eps, float* sfmin, float* tol, aocl_int_t* nsweep, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgsvj1(char* jobv, aocl_int_t* m, aocl_int_t* n, aocl_int_t* n1, scomplex* a, aocl_int_t* lda, scomplex* d, float* sva, aocl_int_t* mv, scomplex* v, aocl_int_t* ldv, float* eps, float* sfmin, float* tol, aocl_int_t* nsweep, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cgtcon(char* norm, aocl_int_t* n, scomplex* dl, scomplex* d, scomplex* du, scomplex* du2, aocl_int_t* ipiv, float* anorm, float* rcond, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cgtrfs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, scomplex* dl, scomplex* d, scomplex* du, scomplex* dlf, scomplex* df, scomplex* duf, scomplex* du2, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgtsv(aocl_int_t* n, aocl_int_t* nrhs, scomplex* dl, scomplex* d, scomplex* du, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cgtsvx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* nrhs, scomplex* dl, scomplex* d, scomplex* du, scomplex* dlf, scomplex* df, scomplex* duf, scomplex* du2, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cgttrf(aocl_int_t* n, scomplex* dl, scomplex* d, scomplex* du, scomplex* du2, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_cgttrs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, scomplex* dl, scomplex* d, scomplex* du, scomplex* du2, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cgtts2(aocl_int_t* itrans, aocl_int_t* n, aocl_int_t* nrhs, scomplex* dl, scomplex* d, scomplex* du, scomplex* du2, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb);
void LAPACK_EXPORT_chb2st_kernels(char* uplo, logical* wantz, aocl_int_t* ttype, aocl_int_t* st, aocl_int_t* ed, aocl_int_t* sweep, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* ib, scomplex* a, aocl_int_t* lda, scomplex* v, scomplex* tau, aocl_int_t* ldvt, scomplex* work);
void LAPACK_EXPORT_chbev_2stage(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_chbev(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_chbevd_2stage(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_chbevd(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_chbevx_2stage(char* jobz, char* range, char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, scomplex* q, aocl_int_t* ldq, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_chbevx(char* jobz, char* range, char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, scomplex* q, aocl_int_t* ldq, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, float* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_chbgst(char* vect, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, scomplex* ab, aocl_int_t* ldab, scomplex* bb, aocl_int_t* ldbb, scomplex* x, aocl_int_t* ldx, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_chbgv(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, scomplex* ab, aocl_int_t* ldab, scomplex* bb, aocl_int_t* ldbb, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_chbgvd(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, scomplex* ab, aocl_int_t* ldab, scomplex* bb, aocl_int_t* ldbb, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_chbgvx(char* jobz, char* range, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, scomplex* ab, aocl_int_t* ldab, scomplex* bb, aocl_int_t* ldbb, scomplex* q, aocl_int_t* ldq, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, float* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_chbtrd(char* vect, char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, float* d, float* e, scomplex* q, aocl_int_t* ldq, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_checon_3(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, float* anorm, float* rcond, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_checon(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, float* anorm, float* rcond, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_checon_rook(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, float* anorm, float* rcond, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cheequb(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* s, float* scond, float* amax, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cheev_2stage(char* jobz, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* w, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cheev(char* jobz, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* w, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cheevd_2stage(char* jobz, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* w, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_cheevd(char* jobz, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* w, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_cheevr_2stage(char* jobz, char* range, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, scomplex* z, aocl_int_t* ldz, aocl_int_t* isuppz, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_cheevr(char* jobz, char* range, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, scomplex* z, aocl_int_t* ldz, aocl_int_t* isuppz, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_cheevx_2stage(char* jobz, char* range, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_cheevx(char* jobz, char* range, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_chegs2(aocl_int_t* itype, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_chegst(aocl_int_t* itype, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_chegv_2stage(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, float* w, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_chegvd(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, float* w, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_chegv(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, float* w, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_chegvx(aocl_int_t* itype, char* jobz, char* range, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_cherfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cherfsx(char* uplo, char* equed, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, float* s, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_chesv_aa_2stage(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chesv_aa(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chesv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chesv_rk(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chesv_rook(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chesvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_chesvxx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, char* equed, float* s, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* rpvgrw, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cheswapr(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* i1, aocl_int_t* i2);
void LAPACK_EXPORT_chetd2(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* d, float* e, scomplex* tau, aocl_int_t* info);
void LAPACK_EXPORT_chetf2(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_chetf2_rk(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_chetf2_rook(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_chetrd_2stage(char* vect, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* d, float* e, scomplex* tau, scomplex* hous2, aocl_int_t* lhous2, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chetrd(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* d, float* e, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chetrd_hb2st(char* stage1, char* vect, char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, float* d, float* e, scomplex* hous, aocl_int_t* lhous, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chetrd_he2hb(char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* a, aocl_int_t* lda, scomplex* ab, aocl_int_t* ldab, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chetrf_aa_2stage(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chetrf_aa(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chetrf(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chetrf_rk(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chetrf_rook(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chetri2(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chetri2x(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* nb, aocl_int_t* info);
void LAPACK_EXPORT_chetri_3(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chetri_3x(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, scomplex* work, aocl_int_t* nb, aocl_int_t* info);
void LAPACK_EXPORT_chetri(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_chetri_rook(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_chetrs2(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_chetrs_3(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_chetrs_aa_2stage(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_chetrs_aa(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_chetrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_chetrs_rook(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_chfrk(char* transr, char* uplo, char* trans, aocl_int_t* n, aocl_int_t* k, float* alpha, scomplex* a, aocl_int_t* lda, float* beta, scomplex* c__);
void LAPACK_EXPORT_chgeqz(char* job, char* compq, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, scomplex* h, aocl_int_t* ldh, scomplex* t, aocl_int_t* ldt, scomplex* alpha, scomplex* beta, scomplex* q, aocl_int_t* ldq, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_chpcon(char* uplo, aocl_int_t* n, scomplex* ap, aocl_int_t* ipiv, float* anorm, float* rcond, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_chpev(char* jobz, char* uplo, aocl_int_t* n, scomplex* ap, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_chpevd(char* jobz, char* uplo, aocl_int_t* n, scomplex* ap, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_chpevx(char* jobz, char* range, char* uplo, aocl_int_t* n, scomplex* ap, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, float* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_chpgst(aocl_int_t* itype, char* uplo, aocl_int_t* n, scomplex* ap, scomplex* bp, aocl_int_t* info);
void LAPACK_EXPORT_chpgvd(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, scomplex* ap, scomplex* bp, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_chpgv(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, scomplex* ap, scomplex* bp, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_chpgvx(aocl_int_t* itype, char* jobz, char* range, char* uplo, aocl_int_t* n, scomplex* ap, scomplex* bp, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, scomplex* z, aocl_int_t* ldz, scomplex* work, float* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_chprfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* ap, scomplex* afp, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_chpsv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* ap, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_chpsvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* ap, scomplex* afp, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_chptrd(char* uplo, aocl_int_t* n, scomplex* ap, float* d, float* e, scomplex* tau, aocl_int_t* info);
void LAPACK_EXPORT_chptrf(char* uplo, aocl_int_t* n, scomplex* ap, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_chptri(char* uplo, aocl_int_t* n, scomplex* ap, aocl_int_t* ipiv, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_chptrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* ap, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_chsein(char* side, char* eigsrc, char* initv, logical* select, aocl_int_t* n, scomplex* h, aocl_int_t* ldh, scomplex* w, scomplex* vl, aocl_int_t* ldvl, scomplex* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, scomplex* work, float* rwork, aocl_int_t* ifaill, aocl_int_t* ifailr, aocl_int_t* info);
void LAPACK_EXPORT_chseqr(char* job, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, scomplex* h, aocl_int_t* ldh, scomplex* w, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_clabrd(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nb, scomplex* a, aocl_int_t* lda, float* d, float* e, scomplex* tauq, scomplex* taup, scomplex* x, aocl_int_t* ldx, scomplex* y, aocl_int_t* ldy);
void LAPACK_EXPORT_clacgv(aocl_int_t* n, scomplex* x, aocl_int_t* incx);
void LAPACK_EXPORT_clacn2(aocl_int_t* n, scomplex* v, scomplex* x, float* est, aocl_int_t* kase, aocl_int_t* isave);
void LAPACK_EXPORT_clacon(aocl_int_t* n, scomplex* v, scomplex* x, float* est, aocl_int_t* kase);
void LAPACK_EXPORT_clacp2(char* uplo, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb);
void LAPACK_EXPORT_clacpy(char* uplo, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb);
void LAPACK_EXPORT_clacrm(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, scomplex* c, aocl_int_t* ldc, float* rwork);
void LAPACK_EXPORT_clacrt(aocl_int_t* n, scomplex* cx, aocl_int_t* incx, scomplex* cy, aocl_int_t* incy, scomplex* c, scomplex* s);
void LAPACK_EXPORT_claed0(aocl_int_t* qsiz, aocl_int_t* n, float* d, float* e, scomplex* q, aocl_int_t* ldq, scomplex* qstore, aocl_int_t* ldqs, float* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_claed7(aocl_int_t* n, aocl_int_t* cutpnt, aocl_int_t* qsiz, aocl_int_t* tlvls, aocl_int_t* curlvl, aocl_int_t* curpbm, float* d, scomplex* q, aocl_int_t* ldq, float* rho, aocl_int_t* indxq, float* qstore, aocl_int_t* qptr, aocl_int_t* prmptr, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, float* givnum, scomplex* work, float* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_claed8(aocl_int_t* k, aocl_int_t* n, aocl_int_t* qsiz, scomplex* q, aocl_int_t* ldq, float* d, float* rho, aocl_int_t* cutpnt, float* z, float* dlamda, scomplex* q2, aocl_int_t* ldq2, float* w, aocl_int_t* indxp, aocl_int_t* indx, aocl_int_t* indxq, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, float* givnum, aocl_int_t* info);
void LAPACK_EXPORT_claein(logical* rightv, logical* noinit, aocl_int_t* n, scomplex* h, aocl_int_t* ldh, scomplex* w, scomplex* v, scomplex* b, aocl_int_t* ldb, float* rwork, float* eps3, float* smlnum, aocl_int_t* info);
void LAPACK_EXPORT_claesy(scomplex* a, scomplex* b, scomplex* c, scomplex* rt1, scomplex* rt2, scomplex* evscal, scomplex* cs1, scomplex* sn1);
void LAPACK_EXPORT_claev2(scomplex* a, scomplex* b, scomplex* c, float* rt1, float* rt2, float* cs1, scomplex* sn1);
void LAPACK_EXPORT_clag2z(aocl_int_t* m, aocl_int_t* n, scomplex* sa, aocl_int_t* ldsa, dcomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_cla_gbamv(aocl_int_t* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, float* alpha, scomplex* ab, aocl_int_t* ldab, scomplex* x, aocl_int_t* incx, float* beta, float* y, aocl_int_t* incy);
void LAPACK_EXPORT_cla_gbrfsx_extended(aocl_int_t* prec_type, aocl_int_t* trans_type, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, scomplex* ab, aocl_int_t* ldab, scomplex* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, logical* colequ, float* c, scomplex* b, aocl_int_t* ldb, scomplex* y, aocl_int_t* ldy, float* berr_out, aocl_int_t* n_norms, float* err_bnds_norm, float* err_bnds_comp, scomplex* res, float* ayb, scomplex* dy, scomplex* y_tail, float* rcond, aocl_int_t* ithresh, float* rthresh, float* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_cla_geamv(aocl_int_t* trans, aocl_int_t* m, aocl_int_t* n, float* alpha, scomplex* a, aocl_int_t* lda, scomplex* x, aocl_int_t* incx, float* beta, float* y, aocl_int_t* incy);
void LAPACK_EXPORT_cla_gerfsx_extended(aocl_int_t* prec_type, aocl_int_t* trans_type, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, logical* colequ, float* c, scomplex* b, aocl_int_t* ldb, scomplex* y, aocl_int_t* ldy, float* berr_out, aocl_int_t* n_norms, float* errs_n, float* errs_c, scomplex* res, float* ayb, scomplex* dy, scomplex* y_tail, float* rcond, aocl_int_t* ithresh, float* rthresh, float* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_clags2(logical* upper, float* a1, scomplex* a2, float* a3, float* b1, scomplex* b2, float* b3, float* csu, scomplex* snu, float* csv, scomplex* snv, float* csq, scomplex* snq);
void LAPACK_EXPORT_clagtm(char* trans, aocl_int_t* n, aocl_int_t* nrhs, float* alpha, scomplex* dl, scomplex* d, scomplex* du, scomplex* x, aocl_int_t* ldx, float* beta, scomplex* b, aocl_int_t* ldb);
void LAPACK_EXPORT_cla_heamv(aocl_int_t* uplo, aocl_int_t* n, float* alpha, scomplex* a, aocl_int_t* lda, scomplex* x, aocl_int_t* incx, float* beta, float* y, aocl_int_t* incy);
void LAPACK_EXPORT_clahef_aa(char* uplo, aocl_int_t* j1, aocl_int_t* m, aocl_int_t* nb, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* h, aocl_int_t* ldh, scomplex* work);
void LAPACK_EXPORT_clahef(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_clahef_rk(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, scomplex* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_clahef_rook(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_cla_herfsx_extended(aocl_int_t* prec_type, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, logical* colequ, float* c, scomplex* b, aocl_int_t* ldb, scomplex* y, aocl_int_t* ldy, float* berr_out, aocl_int_t* n_norms, float* err_bnds_norm, float* err_bnds_comp, scomplex* res, float* ayb, scomplex* dy, scomplex* y_tail, float* rcond, aocl_int_t* ithresh, float* rthresh, float* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_clahqr(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, scomplex* h, aocl_int_t* ldh, scomplex* w, aocl_int_t* iloz, aocl_int_t* ihiz, scomplex* z, aocl_int_t* ldz, aocl_int_t* info);
void LAPACK_EXPORT_clahr2(aocl_int_t* n, aocl_int_t* k, aocl_int_t* nb, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* t, aocl_int_t* ldt, scomplex* y, aocl_int_t* ldy);
void LAPACK_EXPORT_clahrd(aocl_int_t* n, aocl_int_t* k, aocl_int_t* nb, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* t, aocl_int_t* ldt, scomplex* y, aocl_int_t* ldy);
void LAPACK_EXPORT_claic1(aocl_int_t* job, aocl_int_t* j, scomplex* x, float* sest, scomplex* w, scomplex* gamma, float* sestpr, scomplex* s, scomplex* c__);
void LAPACK_EXPORT_cla_lin_berr(aocl_int_t* n, aocl_int_t* nz, aocl_int_t* nrhs, scomplex* res, float* ayb, float* berr);
void LAPACK_EXPORT_clals0(aocl_int_t* icompq, aocl_int_t* nl, aocl_int_t* nr, aocl_int_t* sqre, aocl_int_t* nrhs, scomplex* b, aocl_int_t* ldb, scomplex* bx, aocl_int_t* ldbx, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, aocl_int_t* ldgcol, float* givnum, aocl_int_t* ldgnum, float* poles, float* difl, float* difr, float* z, aocl_int_t* k, float* c, float* s, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_clalsa(aocl_int_t* icompq, aocl_int_t* smlsiz, aocl_int_t* n, aocl_int_t* nrhs, scomplex* b, aocl_int_t* ldb, scomplex* bx, aocl_int_t* ldbx, float* u, aocl_int_t* ldu, float* vt, aocl_int_t* k, float* difl, float* difr, float* z, float* poles, aocl_int_t* givptr, aocl_int_t* givcol, aocl_int_t* ldgcol, aocl_int_t* perm, float* givnum, float* c, float* s, float* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_clalsd(char* uplo, aocl_int_t* smlsiz, aocl_int_t* n, aocl_int_t* nrhs, float* d, float* e, scomplex* b, aocl_int_t* ldb, float* rcond, aocl_int_t* rank, scomplex* work, float* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_clamswlq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* mb, aocl_int_t* nb, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* ldt, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_clamtsqr(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* mb, aocl_int_t* nb, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* ldt, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_clapll(aocl_int_t* n, scomplex* x, aocl_int_t* incx, scomplex* y, aocl_int_t* incy, float* ssmin);
void LAPACK_EXPORT_clapmr(logical* forwrd, aocl_int_t* m, aocl_int_t* n, scomplex* x, aocl_int_t* ldx, aocl_int_t* k);
void LAPACK_EXPORT_clapmt(logical* forwrd, aocl_int_t* m, aocl_int_t* n, scomplex* x, aocl_int_t* ldx, aocl_int_t* k);
void LAPACK_EXPORT_cla_porfsx_extended(aocl_int_t* prec_type, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, logical* colequ, float* c, scomplex* b, aocl_int_t* ldb, scomplex* y, aocl_int_t* ldy, float* berr_out, aocl_int_t* n_norms, float* err_bnds_norm, float* err_bnds_comp, scomplex* res, float* ayb, scomplex* dy, scomplex* y_tail, float* rcond, aocl_int_t* ithresh, float* rthresh, float* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_claqgb(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, scomplex* ab, aocl_int_t* ldab, float* r, float* c, float* rowcnd, float* colcnd, float* amax, char* equed);
void LAPACK_EXPORT_claqge(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, char* equed);
void LAPACK_EXPORT_claqhb(char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, float* s, float* scond, float* amax, char* equed);
void LAPACK_EXPORT_claqhe(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* s, float* scond, float* amax, char* equed);
void LAPACK_EXPORT_claqhp(char* uplo, aocl_int_t* n, scomplex* ap, float* s, float* scond, float* amax, char* equed);
void LAPACK_EXPORT_claqp2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* offset, scomplex* a, aocl_int_t* lda, aocl_int_t* jpvt, scomplex* tau, float* vn1, float* vn2, scomplex* work);
void LAPACK_EXPORT_claqps(aocl_int_t* m, aocl_int_t* n, aocl_int_t* offset, aocl_int_t* nb, aocl_int_t* kb, scomplex* a, aocl_int_t* lda, aocl_int_t* jpvt, scomplex* tau, float* vn1, float* vn2, scomplex* auxv, scomplex* f, aocl_int_t* ldf);
void LAPACK_EXPORT_claqr0(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, scomplex* h, aocl_int_t* ldh, scomplex* w, aocl_int_t* iloz, aocl_int_t* ihiz, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_claqr1(aocl_int_t* n, scomplex* h, aocl_int_t* ldh, scomplex* s1, scomplex* s2, scomplex* v);
void LAPACK_EXPORT_claqr2(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ktop, aocl_int_t* kbot, aocl_int_t* nw, scomplex* h, aocl_int_t* ldh, aocl_int_t* iloz, aocl_int_t* ihiz, scomplex* z, aocl_int_t* ldz, aocl_int_t* ns, aocl_int_t* nd, scomplex* sh, scomplex* v, aocl_int_t* ldv, aocl_int_t* nh, scomplex* t, aocl_int_t* ldt, aocl_int_t* nv, scomplex* wv, aocl_int_t* ldwv, scomplex* work, aocl_int_t* lwork);
void LAPACK_EXPORT_claqr3(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ktop, aocl_int_t* kbot, aocl_int_t* nw, scomplex* h, aocl_int_t* ldh, aocl_int_t* iloz, aocl_int_t* ihiz, scomplex* z, aocl_int_t* ldz, aocl_int_t* ns, aocl_int_t* nd, scomplex* sh, scomplex* v, aocl_int_t* ldv, aocl_int_t* nh, scomplex* t, aocl_int_t* ldt, aocl_int_t* nv, scomplex* wv, aocl_int_t* ldwv, scomplex* work, aocl_int_t* lwork);
void LAPACK_EXPORT_claqr4(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, scomplex* h, aocl_int_t* ldh, scomplex* w, aocl_int_t* iloz, aocl_int_t* ihiz, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_claqr5(logical* wantt, logical* wantz, aocl_int_t* kacc22, aocl_int_t* n, aocl_int_t* ktop, aocl_int_t* kbot, aocl_int_t* nshfts, scomplex* s, scomplex* h, aocl_int_t* ldh, aocl_int_t* iloz, aocl_int_t* ihiz, scomplex* z, aocl_int_t* ldz, scomplex* v, aocl_int_t* ldv, scomplex* u, aocl_int_t* ldu, aocl_int_t* nv, scomplex* wv, aocl_int_t* ldwv, aocl_int_t* nh, scomplex* wh, aocl_int_t* ldwh);
void LAPACK_EXPORT_claqsb(char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, float* s, float* scond, float* amax, char* equed);
void LAPACK_EXPORT_claqsp(char* uplo, aocl_int_t* n, scomplex* ap, float* s, float* scond, float* amax, char* equed);
void LAPACK_EXPORT_claqsy(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* s, float* scond, float* amax, char* equed);
void LAPACK_EXPORT_clar1v(aocl_int_t* n, aocl_int_t* b1, aocl_int_t* bn, float* lambda, float* d, float* l, float* ld, float* lld, float* pivmin, float* gaptol, scomplex* z, logical* wantnc, aocl_int_t* negcnt, float* ztz, float* mingma, aocl_int_t* r, aocl_int_t* isuppz, float* nrminv, float* resid, float* rqcorr, float* work);
void LAPACK_EXPORT_clar2v(aocl_int_t* n, scomplex* x, scomplex* y, scomplex* z, aocl_int_t* incx, float* c, scomplex* s, aocl_int_t* incc);
void LAPACK_EXPORT_clarcm(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* c, aocl_int_t* ldc, float* rwork);
void LAPACK_EXPORT_clarfb(char* side, char* trans, char* direct, char* storev, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* v, aocl_int_t* ldv, scomplex* t, aocl_int_t* ldt, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* ldwork);
void LAPACK_EXPORT_clarf(char* side, aocl_int_t* m, aocl_int_t* n, scomplex* v, aocl_int_t* incv, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work);
void LAPACK_EXPORT_clarfg(aocl_int_t* n, scomplex* alpha, scomplex* x, aocl_int_t* incx, scomplex* tau);
void LAPACK_EXPORT_clarfgp(aocl_int_t* n, scomplex* alpha, scomplex* x, aocl_int_t* incx, scomplex* tau);
void LAPACK_EXPORT_clarft(char* direct, char* storev, aocl_int_t* n, aocl_int_t* k, scomplex* v, aocl_int_t* ldv, scomplex* tau, scomplex* t, aocl_int_t* ldt);
void LAPACK_EXPORT_clarfx(char* side, aocl_int_t* m, aocl_int_t* n, scomplex* v, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work);
void LAPACK_EXPORT_clarfy(char* uplo, aocl_int_t* n, scomplex* v, aocl_int_t* incv, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work);
void LAPACK_EXPORT_clargv(aocl_int_t* n, scomplex* x, aocl_int_t* incx, scomplex* y, aocl_int_t* incy, float* c, aocl_int_t* incc);
void LAPACK_EXPORT_clarnv(aocl_int_t* idist, aocl_int_t* iseed, aocl_int_t* n, scomplex* x);
void LAPACK_EXPORT_clarrv(aocl_int_t* n, float* vl, float* vu, float* d, float* l, float* pivmin, aocl_int_t* isplit, aocl_int_t* m, aocl_int_t* dol, aocl_int_t* dou, float* minrgp, float* rtol1, float* rtol2, float* w, float* werr, float* wgap, aocl_int_t* iblock, aocl_int_t* indexw, float* gers, scomplex* z, aocl_int_t* ldz, aocl_int_t* isuppz, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_clarscl2(aocl_int_t* m, aocl_int_t* n, float* d, scomplex* x, aocl_int_t* ldx);
void LAPACK_EXPORT_clartg(scomplex* f, scomplex* g, float* cs, scomplex* sn, scomplex* r__);
void LAPACK_EXPORT_clartv(aocl_int_t* n, scomplex* x, aocl_int_t* incx, scomplex* y, aocl_int_t* incy, float* c, scomplex* s, aocl_int_t* incc);
void LAPACK_EXPORT_clarzb(char* side, char* trans, char* direct, char* storev, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, scomplex* v, aocl_int_t* ldv, scomplex* t, aocl_int_t* ldt, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* ldwork);
void LAPACK_EXPORT_clarz(char* side, aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, scomplex* v, aocl_int_t* incv, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work);
void LAPACK_EXPORT_clarzt(char* direct, char* storev, aocl_int_t* n, aocl_int_t* k, scomplex* v, aocl_int_t* ldv, scomplex* tau, scomplex* t, aocl_int_t* ldt);
void LAPACK_EXPORT_clascl2(aocl_int_t* m, aocl_int_t* n, float* d, scomplex* x, aocl_int_t* ldx);
void LAPACK_EXPORT_clascl(char* type, aocl_int_t* kl, aocl_int_t* ku, float* cfrom, float* cto, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_claset(char* uplo, aocl_int_t* m, aocl_int_t* n, scomplex* alpha, scomplex* beta, scomplex* a, aocl_int_t* lda);
void LAPACK_EXPORT_clasr(char* side, char* pivot, char* direct, aocl_int_t* m, aocl_int_t* n, float* c, float* s, scomplex* a, aocl_int_t* lda);
void LAPACK_EXPORT_classq(aocl_int_t* n, scomplex* x, aocl_int_t* incx, float* scale, float* sumsq);
void LAPACK_EXPORT_claswlq(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, aocl_int_t* nb, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* ldt, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_claswp(aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* k1, aocl_int_t* k2, aocl_int_t* ipiv, aocl_int_t* incx);
void LAPACK_EXPORT_cla_syamv(aocl_int_t* uplo, aocl_int_t* n, float* alpha, scomplex* a, aocl_int_t* lda, scomplex* x, aocl_int_t* incx, float* beta, float* y, aocl_int_t* incy);
void LAPACK_EXPORT_clasyf_aa(char* uplo, aocl_int_t* j1, aocl_int_t* m, aocl_int_t* nb, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* h, aocl_int_t* ldh, scomplex* work);
void LAPACK_EXPORT_clasyf(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_clasyf_rk(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, scomplex* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_clasyf_rook(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_cla_syrfsx_extended(aocl_int_t* prec_type, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, logical* colequ, float* c, scomplex* b, aocl_int_t* ldb, scomplex* y, aocl_int_t* ldy, float* berr_out, aocl_int_t* n_norms, float* err_bnds_norm, float* err_bnds_comp, scomplex* res, float* ayb, scomplex* dy, scomplex* y_tail, float* rcond, aocl_int_t* ithresh, float* rthresh, float* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_clatbs(char* uplo, char* trans, char* diag, char* normin, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, scomplex* x, float* scale, float* cnorm, aocl_int_t* info);
void LAPACK_EXPORT_clatdf(aocl_int_t* ijob, aocl_int_t* n, scomplex* z, aocl_int_t* ldz, scomplex* rhs, float* rdsum, float* rdscal, aocl_int_t* ipiv, aocl_int_t* jpiv);
void LAPACK_EXPORT_clatps(char* uplo, char* trans, char* diag, char* normin, aocl_int_t* n, scomplex* ap, scomplex* x, float* scale, float* cnorm, aocl_int_t* info);
void LAPACK_EXPORT_clatrd(char* uplo, aocl_int_t* n, aocl_int_t* nb, scomplex* a, aocl_int_t* lda, float* e, scomplex* tau, scomplex* w, aocl_int_t* ldw);
void LAPACK_EXPORT_clatrs(char* uplo, char* trans, char* diag, char* normin, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* x, float* scale, float* cnorm, aocl_int_t* info);
void LAPACK_EXPORT_clatrz(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work);
void LAPACK_EXPORT_clatsqr(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, aocl_int_t* nb, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* ldt, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_clatzm(char* side, aocl_int_t* m, aocl_int_t* n, scomplex* v, aocl_int_t* incv, scomplex* tau, scomplex* c1, scomplex* c2, aocl_int_t* ldc, scomplex* work);
void LAPACK_EXPORT_claunhr_col_getrfnp2(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* d, aocl_int_t* info);
void LAPACK_EXPORT_claunhr_col_getrfnp(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* d, aocl_int_t* info);
void LAPACK_EXPORT_clauu2(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_clauum(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_cla_wwaddw(aocl_int_t* n, scomplex* x, scomplex* y, scomplex* w);
void LAPACK_EXPORT_cpbcon(char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, float* anorm, float* rcond, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cpbequ(char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, float* s, float* scond, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_cpbrfs(char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, scomplex* ab, aocl_int_t* ldab, scomplex* afb, aocl_int_t* ldafb, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cpbstf(char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, aocl_int_t* info);
void LAPACK_EXPORT_cpbsv(char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, scomplex* ab, aocl_int_t* ldab, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cpbsvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, scomplex* ab, aocl_int_t* ldab, scomplex* afb, aocl_int_t* ldafb, char* equed, float* s, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cpbtf2(char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, aocl_int_t* info);
void LAPACK_EXPORT_cpbtrf(char* uplo, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, aocl_int_t* info);
void LAPACK_EXPORT_cpbtrs(char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, scomplex* ab, aocl_int_t* ldab, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cpftrf(char* transr, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* info);
void LAPACK_EXPORT_cpftri(char* transr, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* info);
void LAPACK_EXPORT_cpftrs(char* transr, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cpocon(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* anorm, float* rcond, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cpoequb(aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* s, float* scond, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_cpoequ(aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* s, float* scond, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_cporfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cporfsx(char* uplo, char* equed, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, float* s, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cposv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cposvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, char* equed, float* s, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cposvxx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, char* equed, float* s, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* rpvgrw, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cpotf2(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_cpotrf2(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_cpotrf(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_cpotri(char* uplo, aocl_int_t* n, scomplex* buff_A, aocl_int_t* ldim_A, aocl_int_t* info);
void LAPACK_EXPORT_cpotrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cppcon(char* uplo, aocl_int_t* n, scomplex* ap, float* anorm, float* rcond, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cppequ(char* uplo, aocl_int_t* n, scomplex* ap, float* s, float* scond, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_cpprfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* ap, scomplex* afp, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cppsv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* ap, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cppsvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* ap, scomplex* afp, char* equed, float* s, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cpptrf(char* uplo, aocl_int_t* n, scomplex* ap, aocl_int_t* info);
void LAPACK_EXPORT_cpptri(char* uplo, aocl_int_t* n, scomplex* ap, aocl_int_t* info);
void LAPACK_EXPORT_cpptrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* ap, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cpstf2(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* piv, aocl_int_t* rank, float* tol, float* work, aocl_int_t* info);
void LAPACK_EXPORT_cpstrf(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* piv, aocl_int_t* rank, float* tol, float* work, aocl_int_t* info);
void LAPACK_EXPORT_cptcon(aocl_int_t* n, float* d, scomplex* e, float* anorm, float* rcond, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cpteqr(char* compz, aocl_int_t* n, float* d, float* e, scomplex* z, aocl_int_t* ldz, float* work, aocl_int_t* info);
void LAPACK_EXPORT_cptrfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* d, scomplex* e, float* df, scomplex* ef, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cptsv(aocl_int_t* n, aocl_int_t* nrhs, float* d, scomplex* e, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cptsvx(char* fact, aocl_int_t* n, aocl_int_t* nrhs, float* d, scomplex* e, float* df, scomplex* ef, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cpttrf(aocl_int_t* n, float* d, scomplex* e, aocl_int_t* info);
void LAPACK_EXPORT_cpttrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* d, scomplex* e, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cptts2(aocl_int_t* iuplo, aocl_int_t* n, aocl_int_t* nrhs, float* d, scomplex* e, scomplex* b, aocl_int_t* ldb);
void LAPACK_EXPORT_crot(aocl_int_t* n, scomplex* cx, aocl_int_t* incx, scomplex* cy, aocl_int_t* incy, float* c, scomplex* s);
void LAPACK_EXPORT_cspcon(char* uplo, aocl_int_t* n, scomplex* ap, aocl_int_t* ipiv, float* anorm, float* rcond, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cspmv(char* uplo, aocl_int_t* n, scomplex* alpha, scomplex* ap, scomplex* x, aocl_int_t* incx, scomplex* beta, scomplex* y, aocl_int_t* incy);
void LAPACK_EXPORT_cspr(char* uplo, aocl_int_t* n, scomplex* alpha, scomplex* x, aocl_int_t* incx, scomplex* ap);
void LAPACK_EXPORT_csprfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* ap, scomplex* afp, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_cspsv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* ap, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_cspsvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* ap, scomplex* afp, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_csptrf(char* uplo, aocl_int_t* n, scomplex* ap, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_csptri(char* uplo, aocl_int_t* n, scomplex* ap, aocl_int_t* ipiv, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_csptrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* ap, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_csrscl(aocl_int_t* n, float* sa, scomplex* sx, aocl_int_t* incx);
void LAPACK_EXPORT_cstedc(char* compz, aocl_int_t* n, float* d, float* e, scomplex* z, aocl_int_t* ldz, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_cstegr(char* jobz, char* range, aocl_int_t* n, float* d, float* e, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, scomplex* z, aocl_int_t* ldz, aocl_int_t* isuppz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_cstein(aocl_int_t* n, float* d, float* e, aocl_int_t* m, float* w, aocl_int_t* iblock, aocl_int_t* isplit, scomplex* z, aocl_int_t* ldz, float* work, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_cstemr(char* jobz, char* range, aocl_int_t* n, float* d, float* e, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, aocl_int_t* m, float* w, scomplex* z, aocl_int_t* ldz, aocl_int_t* nzc, aocl_int_t* isuppz, logical* tryrac, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_csteqr(char* compz, aocl_int_t* n, float* d, float* e, scomplex* z, aocl_int_t* ldz, float* work, aocl_int_t* info);
void LAPACK_EXPORT_csteqr(char* jobz, aocl_int_t* n, float* d, float* e, scomplex* z, aocl_int_t* ldz, float* work, aocl_int_t* info);
void LAPACK_EXPORT_csycon_3(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, float* anorm, float* rcond, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_csycon(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, float* anorm, float* rcond, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_csycon_rook(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, float* anorm, float* rcond, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_csyconv(char* uplo, char* way, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_csyconvf(char* uplo, char* way, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_csyconvf_rook(char* uplo, char* way, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_csyequb(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* s, float* scond, float* amax, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_csymv(char* uplo, aocl_int_t* n, scomplex* alpha, scomplex* a, aocl_int_t* lda, scomplex* x, aocl_int_t* incx, scomplex* beta, scomplex* y, aocl_int_t* incy);
void LAPACK_EXPORT_csyr(char* uplo, aocl_int_t* n, scomplex* alpha, scomplex* x, aocl_int_t* incx, scomplex* a, aocl_int_t* lda);
void LAPACK_EXPORT_csyrfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_csyrfsx(char* uplo, char* equed, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, float* s, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_csysv_aa_2stage(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_csysv_aa(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_csysv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_csysv_rk(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_csysv_rook(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_csysvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_csysvxx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, char* equed, float* s, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* rcond, float* rpvgrw, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_csyswapr(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* i1, aocl_int_t* i2);
void LAPACK_EXPORT_csytf2(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_csytf2_rk(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_csytf2_rook(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_csytrf_aa_2stage(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_csytrf_aa(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_csytrf(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_csytrf_rk(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_csytrf_rook(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_csytri2(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_csytri2x(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* nb, aocl_int_t* info);
void LAPACK_EXPORT_csytri_3(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_csytri_3x(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, scomplex* work, aocl_int_t* nb, aocl_int_t* info);
void LAPACK_EXPORT_csytri(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_csytri_rook(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_csytrs2(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_csytrs_3(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* e, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_csytrs_aa_2stage(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_csytrs_aa(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_csytrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_csytrs_rook(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_ctbcon(char* norm, char* uplo, char* diag, aocl_int_t* n, aocl_int_t* kd, scomplex* ab, aocl_int_t* ldab, float* rcond, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ctbrfs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, scomplex* ab, aocl_int_t* ldab, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ctbtrs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, scomplex* ab, aocl_int_t* ldab, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_ctfsm(char* transr, char* side, char* uplo, char* trans, char* diag, aocl_int_t* m, aocl_int_t* n, scomplex* alpha, scomplex* a, scomplex* b, aocl_int_t* ldb);
void LAPACK_EXPORT_ctftri(char* transr, char* uplo, char* diag, aocl_int_t* n, scomplex* a, aocl_int_t* info);
void LAPACK_EXPORT_ctfttp(char* transr, char* uplo, aocl_int_t* n, scomplex* arf, scomplex* ap, aocl_int_t* info);
void LAPACK_EXPORT_ctfttr(char* transr, char* uplo, aocl_int_t* n, scomplex* arf, scomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_ctgevc(char* side, char* howmny, logical* select, aocl_int_t* n, scomplex* s, aocl_int_t* lds, scomplex* p, aocl_int_t* ldp, scomplex* vl, aocl_int_t* ldvl, scomplex* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ctgex2(logical* wantq, logical* wantz, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* q, aocl_int_t* ldq, scomplex* z, aocl_int_t* ldz, aocl_int_t* j1, aocl_int_t* info);
void LAPACK_EXPORT_ctgexc(logical* wantq, logical* wantz, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* q, aocl_int_t* ldq, scomplex* z, aocl_int_t* ldz, aocl_int_t* ifst, aocl_int_t* ilst, aocl_int_t* info);
void LAPACK_EXPORT_ctgsen(aocl_int_t* ijob, logical* wantq, logical* wantz, logical* select, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* alpha, scomplex* beta, scomplex* q, aocl_int_t* ldq, scomplex* z, aocl_int_t* ldz, aocl_int_t* m, float* pl, float* pr, float* dif, scomplex* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_ctgsja(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, float* tola, float* tolb, float* alpha, float* beta, scomplex* u, aocl_int_t* ldu, scomplex* v, aocl_int_t* ldv, scomplex* q, aocl_int_t* ldq, scomplex* work, aocl_int_t* ncycle, aocl_int_t* info);
void LAPACK_EXPORT_ctgsna(char* job, char* howmny, logical* select, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* vl, aocl_int_t* ldvl, scomplex* vr, aocl_int_t* ldvr, float* s, float* dif, aocl_int_t* mm, aocl_int_t* m, scomplex* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_ctgsy2(char* trans, aocl_int_t* ijob, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* c, aocl_int_t* ldc, scomplex* d, aocl_int_t* ldd, scomplex* e, aocl_int_t* lde, scomplex* f, aocl_int_t* ldf, float* scale, float* rdsum, float* rdscal, aocl_int_t* info);
void LAPACK_EXPORT_ctgsyl(char* trans, aocl_int_t* ijob, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* c, aocl_int_t* ldc, scomplex* d, aocl_int_t* ldd, scomplex* e, aocl_int_t* lde, scomplex* f, aocl_int_t* ldf, float* scale, float* dif, scomplex* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_ctpcon(char* norm, char* uplo, char* diag, aocl_int_t* n, scomplex* ap, float* rcond, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ctplqt2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_ctplqt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, aocl_int_t* mb, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* t, aocl_int_t* ldt, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_ctpmlqt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, aocl_int_t* mb, scomplex* v, aocl_int_t* ldv, scomplex* t, aocl_int_t* ldt, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_ctpmqrt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, aocl_int_t* nb, scomplex* v, aocl_int_t* ldv, scomplex* t, aocl_int_t* ldt, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_ctpqrt2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_ctpqrt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, aocl_int_t* nb, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* t, aocl_int_t* ldt, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_ctprfb(char* side, char* trans, char* direct, char* storev, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, scomplex* v, aocl_int_t* ldv, scomplex* t, aocl_int_t* ldt, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* work, aocl_int_t* ldwork);
void LAPACK_EXPORT_ctprfs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, scomplex* ap, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ctptri(char* uplo, char* diag, aocl_int_t* n, scomplex* ap, aocl_int_t* info);
void LAPACK_EXPORT_ctptrs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, scomplex* ap, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_ctpttf(char* transr, char* uplo, aocl_int_t* n, scomplex* ap, scomplex* arf, aocl_int_t* info);
void LAPACK_EXPORT_ctpttr(char* uplo, aocl_int_t* n, scomplex* ap, scomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_ctrcon(char* norm, char* uplo, char* diag, aocl_int_t* n, scomplex* a, aocl_int_t* lda, float* rcond, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ctrevc3(char* side, char* howmny, logical* select, aocl_int_t* n, scomplex* t, aocl_int_t* ldt, scomplex* vl, aocl_int_t* ldvl, scomplex* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* info);
void LAPACK_EXPORT_ctrevc(char* side, char* howmny, logical* select, aocl_int_t* n, scomplex* t, aocl_int_t* ldt, scomplex* vl, aocl_int_t* ldvl, scomplex* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ctrexc(char* compq, aocl_int_t* n, scomplex* t, aocl_int_t* ldt, scomplex* q, aocl_int_t* ldq, aocl_int_t* ifst, aocl_int_t* ilst, aocl_int_t* info);
void LAPACK_EXPORT_ctrrfs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* x, aocl_int_t* ldx, float* ferr, float* berr, scomplex* work, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ctrsen(char* job, char* compq, logical* select, aocl_int_t* n, scomplex* t, aocl_int_t* ldt, scomplex* q, aocl_int_t* ldq, scomplex* w, aocl_int_t* m, float* s, float* sep, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ctrsna(char* job, char* howmny, logical* select, aocl_int_t* n, scomplex* t, aocl_int_t* ldt, scomplex* vl, aocl_int_t* ldvl, scomplex* vr, aocl_int_t* ldvr, float* s, float* sep, aocl_int_t* mm, aocl_int_t* m, scomplex* work, aocl_int_t* ldwork, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ctrsyl(char* trana, char* tranb, aocl_int_t* isgn, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, scomplex* c, aocl_int_t* ldc, float* scale, aocl_int_t* info);
void LAPACK_EXPORT_ctrti2(char* uplo, char* diag, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_ctrtri(char* uplo, char* diag, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_ctrtrs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, scomplex* a, aocl_int_t* lda, scomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_ctrttf(char* transr, char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* arf, aocl_int_t* info);
void LAPACK_EXPORT_ctrttp(char* uplo, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* ap, aocl_int_t* info);
void LAPACK_EXPORT_ctzrqf(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tau, aocl_int_t* info);
void LAPACK_EXPORT_ctzrzf(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunbdb1(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, scomplex* x11, aocl_int_t* ldx11, scomplex* x21, aocl_int_t* ldx21, float* theta, float* phi, scomplex* taup1, scomplex* taup2, scomplex* tauq1, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunbdb2(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, scomplex* x11, aocl_int_t* ldx11, scomplex* x21, aocl_int_t* ldx21, float* theta, float* phi, scomplex* taup1, scomplex* taup2, scomplex* tauq1, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunbdb3(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, scomplex* x11, aocl_int_t* ldx11, scomplex* x21, aocl_int_t* ldx21, float* theta, float* phi, scomplex* taup1, scomplex* taup2, scomplex* tauq1, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunbdb4(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, scomplex* x11, aocl_int_t* ldx11, scomplex* x21, aocl_int_t* ldx21, float* theta, float* phi, scomplex* taup1, scomplex* taup2, scomplex* tauq1, scomplex* phantom, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunbdb5(aocl_int_t* m1, aocl_int_t* m2, aocl_int_t* n, scomplex* x1, aocl_int_t* incx1, scomplex* x2, aocl_int_t* incx2, scomplex* q1, aocl_int_t* ldq1, scomplex* q2, aocl_int_t* ldq2, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunbdb6(aocl_int_t* m1, aocl_int_t* m2, aocl_int_t* n, scomplex* x1, aocl_int_t* incx1, scomplex* x2, aocl_int_t* incx2, scomplex* q1, aocl_int_t* ldq1, scomplex* q2, aocl_int_t* ldq2, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunbdb(char* trans, char* signs, aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, scomplex* x11, aocl_int_t* ldx11, scomplex* x12, aocl_int_t* ldx12, scomplex* x21, aocl_int_t* ldx21, scomplex* x22, aocl_int_t* ldx22, float* theta, float* phi, scomplex* taup1, scomplex* taup2, scomplex* tauq1, scomplex* tauq2, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cuncsd2by1(char* jobu1, char* jobu2, char* jobv1t, aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, scomplex* x11, aocl_int_t* ldx11, scomplex* x21, aocl_int_t* ldx21, float* theta, scomplex* u1, aocl_int_t* ldu1, scomplex* u2, aocl_int_t* ldu2, scomplex* v1t, aocl_int_t* ldv1t, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_cuncsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, char* signs, aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, scomplex* x11, aocl_int_t* ldx11, scomplex* x12, aocl_int_t* ldx12, scomplex* x21, aocl_int_t* ldx21, scomplex* x22, aocl_int_t* ldx22, float* theta, scomplex* u1, aocl_int_t* ldu1, scomplex* u2, aocl_int_t* ldu2, scomplex* v1t, aocl_int_t* ldv1t, scomplex* v2t, aocl_int_t* ldv2t, scomplex* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_cung2l(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cung2r(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cungbr(char* vect, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunghr(aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cungl2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cunglq(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cungql(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cungqr(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cungr2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cungrq(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cungtr(char* uplo, aocl_int_t* m, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cungtsqr(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, aocl_int_t* nb, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* ldt, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunhr_col(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nb, scomplex* a, aocl_int_t* lda, scomplex* t, aocl_int_t* ldt, scomplex* d, aocl_int_t* info);
void LAPACK_EXPORT_cunm22(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* n1, aocl_int_t* n2, scomplex* q, aocl_int_t* ldq, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunm2l(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cunm2r(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cunmbr(char* vect, char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunmhr(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunml2(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cunmlq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunmql(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunmqr(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunmr2(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cunmr3(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cunmrq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunmrz(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cunmtr(char* side, char* uplo, char* trans, aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_cupgtr(char* uplo, aocl_int_t* n, scomplex* ap, scomplex* tau, scomplex* q, aocl_int_t* ldq, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_cupmtr(char* side, char* uplo, char* trans, aocl_int_t* m, aocl_int_t* n, scomplex* ap, scomplex* tau, scomplex* c, aocl_int_t* ldc, scomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_dbbcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, double* theta, double* phi, double* u1, aocl_int_t* ldu1, double* u2, aocl_int_t* ldu2, double* v1t, aocl_int_t* ldv1t, double* v2t, aocl_int_t* ldv2t, double* b11d, double* b11e, double* b12d, double* b12e, double* b21d, double* b21e, double* b22d, double* b22e, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dbdsdc(char* uplo, char* compq, aocl_int_t* n, double* d, double* e, double* u, aocl_int_t* ldu, double* vt, aocl_int_t* ldvt, double* q, aocl_int_t* iq, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dbdsqr(char* uplo, aocl_int_t* n, aocl_int_t* ncvt, aocl_int_t* nru, aocl_int_t* ncc, double* d, double* e, double* vt, aocl_int_t* ldvt, double* u, aocl_int_t* ldu, double* c, aocl_int_t* ldc, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_dbdsvdx(char* uplo, char* jobz, char* range, aocl_int_t* n, double* d, double* e, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, aocl_int_t* ns, double* s, double* z, aocl_int_t* ldz, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dcombssq(double* v1, double* v2);
void LAPACK_EXPORT_ddisna(char* job, aocl_int_t* m, aocl_int_t* n, double* d, double* sep, aocl_int_t* info);
void LAPACK_EXPORT_dgbbrd(char* vect, aocl_int_t* m, aocl_int_t* n, aocl_int_t* ncc, aocl_int_t* kl, aocl_int_t* ku, double* ab, aocl_int_t* ldab, double* d, double* e, double* q, aocl_int_t* ldq, double* pt, aocl_int_t* ldpt, double* c, aocl_int_t* ldc, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgbcon(char* norm, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, double* ab, aocl_int_t* ldab, aocl_int_t* ipiv, double* anorm, double* rcond, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgbequb(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, double* ab, aocl_int_t* ldab, double* r, double* c, double* rowcnd, double* colcnd, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_dgbequ(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, double* ab, aocl_int_t* ldab, double* r, double* c, double* rowcnd, double* colcnd, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_dgbrfs(char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, double* ab, aocl_int_t* ldab, double* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgbrfsx(char* trans, char* equed, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, double* ab, aocl_int_t* ldab, double* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, double* r, double* c, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgbsv(aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, double* ab, aocl_int_t* ldab, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dgbsvx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, double* ab, aocl_int_t* ldab, double* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, char* equed, double* r, double* c, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgbsvxx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, double* ab, aocl_int_t* ldab, double* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, char* equed, double* r, double* c, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* rpvgrw, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgbtf2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, double* ab, aocl_int_t* ldab, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_dgbtrf(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, double* ab, aocl_int_t* ldab, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_dgbtrs(char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, double* ab, aocl_int_t* ldab, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dgebak(char* job, char* side, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, double* scale, aocl_int_t* m, double* v, aocl_int_t* ldv, aocl_int_t* info);
void LAPACK_EXPORT_dgebal(char* job, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ilo, aocl_int_t* ihi, double* scale, aocl_int_t* info);
void LAPACK_EXPORT_dgebd2(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* d, double* e, double* tauq, double* taup, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgebrd(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* d, double* e, double* tauq, double* taup, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgecon(char* norm, aocl_int_t* n, double* a, aocl_int_t* lda, double* anorm, double* rcond, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgeequb(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_dgeequ(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_dgees(char* jobvs, char* sort, L_fp select, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* sdim, double* wr, double* wi, double* vs, aocl_int_t* ldvs, double* work, aocl_int_t* lwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_dgeesx(char* jobvs, char* sort, L_fp select, char* sense, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* sdim, double* wr, double* wi, double* vs, aocl_int_t* ldvs, double* rconde, double* rcondv, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_dgeev(char* jobvl, char* jobvr, aocl_int_t* n, double* a, aocl_int_t* lda, double* wr, double* wi, double* vl, aocl_int_t* ldvl, double* vr, aocl_int_t* ldvr, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgeevx(char* balanc, char* jobvl, char* jobvr, char* sense, aocl_int_t* n, double* a, aocl_int_t* lda, double* wr, double* wi, double* vl, aocl_int_t* ldvl, double* vr, aocl_int_t* ldvr, aocl_int_t* ilo, aocl_int_t* ihi, double* scale, double* abnrm, double* rconde, double* rcondv, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgegs(char* jobvsl, char* jobvsr, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* alphar, double* alphai, double* beta, double* vsl, aocl_int_t* ldvsl, double* vsr, aocl_int_t* ldvsr, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgegv(char* jobvl, char* jobvr, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* alphar, double* alphai, double* beta, double* vl, aocl_int_t* ldvl, double* vr, aocl_int_t* ldvr, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgehd2(aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgehrd(aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgejsv(char* joba, char* jobu, char* jobv, char* jobr, char* jobt, char* jobp, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* sva, double* u, aocl_int_t* ldu, double* v, aocl_int_t* ldv, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgelq2(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgelqf(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgelq(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* t, aocl_int_t* tsize, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgelqt3(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_dgelqt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, double* a, aocl_int_t* lda, double* t, aocl_int_t* ldt, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgels(char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgelsd(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* s, double* rcond, aocl_int_t* rank, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgelss(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* s, double* rcond, aocl_int_t* rank, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgelsx(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, aocl_int_t* jpvt, double* rcond, aocl_int_t* rank, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgelsy(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, aocl_int_t* jpvt, double* rcond, aocl_int_t* rank, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgemlq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* t, aocl_int_t* tsize, double* c, aocl_int_t* ldc, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgemlqt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* mb, double* v, aocl_int_t* ldv, double* t, aocl_int_t* ldt, double* c, aocl_int_t* ldc, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgemqr(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* t, aocl_int_t* tsize, double* c, aocl_int_t* ldc, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgemqrt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* nb, double* v, aocl_int_t* ldv, double* t, aocl_int_t* ldt, double* c, aocl_int_t* ldc, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgeql2(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgeqlf(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgeqp3(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* jpvt, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgeqpf(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* jpvt, double* tau, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgeqr2(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgeqr2p(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgeqrf(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgeqrfp(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgeqr(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* t, aocl_int_t* tsize, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgeqrt2(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_dgeqrt3(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_dgeqrt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nb, double* a, aocl_int_t* lda, double* t, aocl_int_t* ldt, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgerfs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgerfsx(char* trans, char* equed, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, aocl_int_t* ipiv, double* r, double* c, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgerq2(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgerqf(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgesc2(aocl_int_t* n, double* a, aocl_int_t* lda, double* rhs, aocl_int_t* ipiv, aocl_int_t* jpiv, double* scale);
void LAPACK_EXPORT_dgesdd(char* jobz, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* s, double* u, aocl_int_t* ldu, double* vt, aocl_int_t* ldvt, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgesvd(char* jobu, char* jobv, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* s, double* u, aocl_int_t* ldu, double* vt, aocl_int_t* ldvt, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgesvdq(char* joba, char* jobp, char* jobr, char* jobu, char* jobv, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* s, double* u, aocl_int_t* ldu, double* v, aocl_int_t* ldv, aocl_int_t* numrank, aocl_int_t* iwork, aocl_int_t* liwork, double* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* info);
void LAPACK_EXPORT_dgesvdx(char* jobu, char* jobvt, char* range, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, aocl_int_t* ns, double* s, double* u, aocl_int_t* ldu, double* vt, aocl_int_t* ldvt, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgesv(aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dgesvj(char* joba, char* jobu, char* jobv, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* sva, aocl_int_t* mv, double* v, aocl_int_t* ldv, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgesvx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, aocl_int_t* ipiv, char* equed, double* r, double* c, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgesvxx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, aocl_int_t* ipiv, char* equed, double* r, double* c, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* rpvgrw, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgetc2(aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* jpiv, aocl_int_t* info);
void LAPACK_EXPORT_dgetf2(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_dgetrf2(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_dgetrf(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_dgetri(aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgetrs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dgetsls(char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dggbak(char* job, char* side, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, double* lscale, double* rscale, aocl_int_t* m, double* v, aocl_int_t* ldv, aocl_int_t* info);
void LAPACK_EXPORT_dggbal(char* job, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, aocl_int_t* ilo, aocl_int_t* ihi, double* lscale, double* rscale, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgges3(char* jobvsl, char* jobvsr, char* sort, L_fp selctg, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, aocl_int_t* sdim, double* alphar, double* alphai, double* beta, double* vsl, aocl_int_t* ldvsl, double* vsr, aocl_int_t* ldvsr, double* work, aocl_int_t* lwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_dgges(char* jobvsl, char* jobvsr, char* sort, L_fp selctg, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, aocl_int_t* sdim, double* alphar, double* alphai, double* beta, double* vsl, aocl_int_t* ldvsl, double* vsr, aocl_int_t* ldvsr, double* work, aocl_int_t* lwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_dggesx(char* jobvsl, char* jobvsr, char* sort, L_fp selctg, char* sense, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, aocl_int_t* sdim, double* alphar, double* alphai, double* beta, double* vsl, aocl_int_t* ldvsl, double* vsr, aocl_int_t* ldvsr, double* rconde, double* rcondv, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_dggev3(char* jobvl, char* jobvr, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* alphar, double* alphai, double* beta, double* vl, aocl_int_t* ldvl, double* vr, aocl_int_t* ldvr, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dggev(char* jobvl, char* jobvr, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* alphar, double* alphai, double* beta, double* vl, aocl_int_t* ldvl, double* vr, aocl_int_t* ldvr, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dggevx(char* balanc, char* jobvl, char* jobvr, char* sense, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* alphar, double* alphai, double* beta, double* vl, aocl_int_t* ldvl, double* vr, aocl_int_t* ldvr, aocl_int_t* ilo, aocl_int_t* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, double* rconde, double* rcondv, double* work, aocl_int_t* lwork, aocl_int_t* iwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_dggglm(aocl_int_t* n, aocl_int_t* m, aocl_int_t* p, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* d, double* x, double* y, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgghd3(char* compq, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* q, aocl_int_t* ldq, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgghrd(char* compq, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* q, aocl_int_t* ldq, double* z, aocl_int_t* ldz, aocl_int_t* info);
void LAPACK_EXPORT_dgglse(aocl_int_t* m, aocl_int_t* n, aocl_int_t* p, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* c, double* d, double* x, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dggqrf(aocl_int_t* n, aocl_int_t* m, aocl_int_t* p, double* a, aocl_int_t* lda, double* taua, double* b, aocl_int_t* ldb, double* taub, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dggrqf(aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, double* a, aocl_int_t* lda, double* taua, double* b, aocl_int_t* ldb, double* taub, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dggsvd3(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* n, aocl_int_t* p, aocl_int_t* k, aocl_int_t* l, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* alpha, double* beta, double* u, aocl_int_t* ldu, double* v, aocl_int_t* ldv, double* q, aocl_int_t* ldq, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dggsvd(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* n, aocl_int_t* p, aocl_int_t* k, aocl_int_t* l, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* alpha, double* beta, double* u, aocl_int_t* ldu, double* v, aocl_int_t* ldv, double* q, aocl_int_t* ldq, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dggsvp3(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* tola, double* tolb, aocl_int_t* k, aocl_int_t* l, double* u, aocl_int_t* ldu, double* v, aocl_int_t* ldv, double* q, aocl_int_t* ldq, aocl_int_t* iwork, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dggsvp(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* tola, double* tolb, aocl_int_t* k, aocl_int_t* l, double* u, aocl_int_t* ldu, double* v, aocl_int_t* ldv, double* q, aocl_int_t* ldq, aocl_int_t* iwork, double* tau, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dgsvj0(char* jobv, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* d, double* sva, aocl_int_t* mv, double* v, aocl_int_t* ldv, double* eps, double* sfmin, double* tol, aocl_int_t* nsweep, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgsvj1(char* jobv, aocl_int_t* m, aocl_int_t* n, aocl_int_t* n1, double* a, aocl_int_t* lda, double* d, double* sva, aocl_int_t* mv, double* v, aocl_int_t* ldv, double* eps, double* sfmin, double* tol, aocl_int_t* nsweep, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dgtcon(char* norm, aocl_int_t* n, double* dl, double* d, double* du, double* du2, aocl_int_t* ipiv, double* anorm, double* rcond, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgtrfs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, double* dl, double* d, double* du, double* dlf, double* df, double* duf, double* du2, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgtsv(aocl_int_t* n, aocl_int_t* nrhs, double* dl, double* d, double* du, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dgtsvx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* nrhs, double* dl, double* d, double* du, double* dlf, double* df, double* duf, double* du2, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dgttrf(aocl_int_t* n, double* dl, double* d, double* du, double* du2, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_dgttrs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, double* dl, double* d, double* du, double* du2, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dgtts2(aocl_int_t* itrans, aocl_int_t* n, aocl_int_t* nrhs, double* dl, double* d, double* du, double* du2, aocl_int_t* ipiv, double* b, aocl_int_t* ldb);
void LAPACK_EXPORT_dhgeqz(char* job, char* compq, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, double* h, aocl_int_t* ldh, double* t, aocl_int_t* ldt, double* alphar, double* alphai, double* beta, double* q, aocl_int_t* ldq, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dhsein(char* side, char* eigsrc, char* initv, logical* select, aocl_int_t* n, double* h, aocl_int_t* ldh, double* wr, double* wi, double* vl, aocl_int_t* ldvl, double* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, double* work, aocl_int_t* ifaill, aocl_int_t* ifailr, aocl_int_t* info);
void LAPACK_EXPORT_dhseqr(char* job, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, double* h, aocl_int_t* ldh, double* wr, double* wi, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dlabad(double* small_, double* large);
void LAPACK_EXPORT_dlabrd(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nb, double* a, aocl_int_t* lda, double* d, double* e, double* tauq, double* taup, double* x, aocl_int_t* ldx, double* y, aocl_int_t* ldy);
void LAPACK_EXPORT_dlacn2(aocl_int_t* n, double* v, double* x, aocl_int_t* isgn, double* est, aocl_int_t* kase, aocl_int_t* isave);
void LAPACK_EXPORT_dlacon(aocl_int_t* n, double* v, double* x, aocl_int_t* isgn, double* est, aocl_int_t* kase);
void LAPACK_EXPORT_dlacpy(char* uplo, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb);
void LAPACK_EXPORT_dladiv(double* a, double* b, double* c, double* d, double* p, double* q);
void LAPACK_EXPORT_dlae2(double* a, double* b, double* c, double* rt1, double* rt2);
void LAPACK_EXPORT_dlaebz(aocl_int_t* ijob, aocl_int_t* nitmax, aocl_int_t* n, aocl_int_t* mmax, aocl_int_t* minp, aocl_int_t* nbmin, double* abstol, double* reltol, double* pivmin, double* d, double* e, double* e2, aocl_int_t* nval, double* ab, double* c, aocl_int_t* mout, aocl_int_t* nab, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dlaed0(aocl_int_t* icompq, aocl_int_t* qsiz, aocl_int_t* n, double* d, double* e, double* q, aocl_int_t* ldq, double* qstore, aocl_int_t* ldqs, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dlaed1(aocl_int_t* n, double* d, double* q, aocl_int_t* ldq, aocl_int_t* indxq, double* rho, aocl_int_t* cutpnt, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dlaed2(aocl_int_t* k, aocl_int_t* n, aocl_int_t* n1, double* d, double* q, aocl_int_t* ldq, aocl_int_t* indxq, double* rho, double* z, double* dlamda, double* w, double* q2, aocl_int_t* indx, aocl_int_t* indxc, aocl_int_t* indxp, aocl_int_t* coltyp, aocl_int_t* info);
void LAPACK_EXPORT_dlaed3(aocl_int_t* k, aocl_int_t* n, aocl_int_t* n1, double* d, double* q, aocl_int_t* ldq, double* rho, double* dlamda, double* q2, aocl_int_t* indx, aocl_int_t* ctot, double* w, double* s, aocl_int_t* info);
void LAPACK_EXPORT_dlaed4(aocl_int_t* n, aocl_int_t* i, double* d, double* z, double* delta, double* rho, double* dlam, aocl_int_t* info);
void LAPACK_EXPORT_dlaed5(aocl_int_t* i, double* d, double* z, double* delta, double* rho, double* dlam);
void LAPACK_EXPORT_dlaed6(aocl_int_t* kniter, logical* orgati, double* rho, double* d, double* z, double* finit, double* tau, aocl_int_t* info);
void LAPACK_EXPORT_dlaed7(aocl_int_t* icompq, aocl_int_t* n, aocl_int_t* qsiz, aocl_int_t* tlvls, aocl_int_t* curlvl, aocl_int_t* curpbm, double* d, double* q, aocl_int_t* ldq, aocl_int_t* indxq, double* rho, aocl_int_t* cutpnt, double* qstore, aocl_int_t* qptr, aocl_int_t* prmptr, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, double* givnum, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dlaed8(aocl_int_t* icompq, aocl_int_t* k, aocl_int_t* n, aocl_int_t* qsiz, double* d, double* q, aocl_int_t* ldq, aocl_int_t* indxq, double* rho, aocl_int_t* cutpnt, double* z, double* dlamda, double* q2, aocl_int_t* ldq2, double* w, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, double* givnum, aocl_int_t* indxp, aocl_int_t* indx, aocl_int_t* info);
void LAPACK_EXPORT_dlaed9(aocl_int_t* k, aocl_int_t* kstart, aocl_int_t* kstop, aocl_int_t* n, double* d, double* q, aocl_int_t* ldq, double* rho, double* dlamda, double* w, double* s, aocl_int_t* lds, aocl_int_t* info);
void LAPACK_EXPORT_dlaeda(aocl_int_t* n, aocl_int_t* tlvls, aocl_int_t* curlvl, aocl_int_t* curpbm, aocl_int_t* prmptr, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, double* givnum, double* q, aocl_int_t* qptr, double* z, double* ztemp, aocl_int_t* info);
void LAPACK_EXPORT_dlaein(logical* rightv, logical* noinit, aocl_int_t* n, double* h, aocl_int_t* ldh, double* wr, double* wi, double* vr, double* vi, double* b, aocl_int_t* ldb, double* work, double* eps3, double* smlnum, double* bignum, aocl_int_t* info);
void LAPACK_EXPORT_dlaev2(double* a, double* b, double* c, double* rt1, double* rt2, double* cs1, double* sn1);
void LAPACK_EXPORT_dlaexc(logical* wantq, aocl_int_t* n, double* t, aocl_int_t* ldt, double* q, aocl_int_t* ldq, aocl_int_t* j1, aocl_int_t* n1, aocl_int_t* n2, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dlag2(double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* safmin, double* scale1, double* scale2, double* wr1, double* wr2, double* wi);
void LAPACK_EXPORT_dlag2s(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, float* sa, aocl_int_t* ldsa, aocl_int_t* info);
void LAPACK_EXPORT_dla_gbamv(aocl_int_t* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, double* alpha, double* ab, aocl_int_t* ldab, double* x, aocl_int_t* incx, double* beta, double* y, aocl_int_t* incy);
void LAPACK_EXPORT_dla_gbrfsx_extended(aocl_int_t* prec_type, aocl_int_t* trans_type, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, double* ab, aocl_int_t* ldab, double* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, logical* colequ, double* c, double* b, aocl_int_t* ldb, double* y, aocl_int_t* ldy, double* berr_out, aocl_int_t* n_norms, double* err_bnds_norm, double* err_bnds_comp, double* res, double* ayb, double* dy, double* y_tail, double* rcond, aocl_int_t* ithresh, double* rthresh, double* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_dla_geamv(aocl_int_t* trans, aocl_int_t* m, aocl_int_t* n, double* alpha, double* a, aocl_int_t* lda, double* x, aocl_int_t* incx, double* beta, double* y, aocl_int_t* incy);
void LAPACK_EXPORT_dla_gerfsx_extended(aocl_int_t* prec_type, aocl_int_t* trans_type, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, aocl_int_t* ipiv, logical* colequ, double* c, double* b, aocl_int_t* ldb, double* y, aocl_int_t* ldy, double* berr_out, aocl_int_t* n_norms, double* errs_n, double* errs_c, double* res, double* ayb, double* dy, double* y_tail, double* rcond, aocl_int_t* ithresh, double* rthresh, double* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_dlags2(logical* upper, double* a1, double* a2, double* a3, double* b1, double* b2, double* b3, double* csu, double* snu, double* csv, double* snv, double* csq, double* snq);
void LAPACK_EXPORT_dlagtf(aocl_int_t* n, double* a, double* lambda, double* b, double* c, double* tol, double* d, aocl_int_t* in, aocl_int_t* info);
void LAPACK_EXPORT_dlagtm(char* trans, aocl_int_t* n, aocl_int_t* nrhs, double* alpha, double* dl, double* d, double* du, double* x, aocl_int_t* ldx, double* beta, double* b, aocl_int_t* ldb);
void LAPACK_EXPORT_dlagts(aocl_int_t* job, aocl_int_t* n, double* a, double* b, double* c, double* d, aocl_int_t* in, double* y, double* tol, aocl_int_t* info);
void LAPACK_EXPORT_dlagv2(double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* alphar, double* alphai, double* beta, double* csl, double* snl, double* csr, double* snr);
void LAPACK_EXPORT_dlahqr(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, double* h, aocl_int_t* ldh, double* wr, double* wi, aocl_int_t* iloz, aocl_int_t* ihiz, double* z, aocl_int_t* ldz, aocl_int_t* info);
void LAPACK_EXPORT_dlahr2(aocl_int_t* n, aocl_int_t* k, aocl_int_t* nb, double* a, aocl_int_t* lda, double* tau, double* t, aocl_int_t* ldt, double* y, aocl_int_t* ldy);
void LAPACK_EXPORT_dlahrd(aocl_int_t* n, aocl_int_t* k, aocl_int_t* nb, double* a, aocl_int_t* lda, double* tau, double* t, aocl_int_t* ldt, double* y, aocl_int_t* ldy);
void LAPACK_EXPORT_dlaic1(aocl_int_t* job, aocl_int_t* j, double* x, double* sest, double* w, double* gamma, double* sestpr, double* s, double* c__);
void LAPACK_EXPORT_dla_lin_berr(aocl_int_t* n, aocl_int_t* nz, aocl_int_t* nrhs, double* res, double* ayb, double* berr);
void LAPACK_EXPORT_dlaln2(logical* ltrans, aocl_int_t* na, aocl_int_t* nw, double* smin, double* ca, double* a, aocl_int_t* lda, double* d1, double* d2, double* b, aocl_int_t* ldb, double* wr, double* wi, double* x, aocl_int_t* ldx, double* scale, double* xnorm, aocl_int_t* info);
void LAPACK_EXPORT_dlals0(aocl_int_t* icompq, aocl_int_t* nl, aocl_int_t* nr, aocl_int_t* sqre, aocl_int_t* nrhs, double* b, aocl_int_t* ldb, double* bx, aocl_int_t* ldbx, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, aocl_int_t* ldgcol, double* givnum, aocl_int_t* ldgnum, double* poles, double* difl, double* difr, double* z, aocl_int_t* k, double* c, double* s, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dlalsa(aocl_int_t* icompq, aocl_int_t* smlsiz, aocl_int_t* n, aocl_int_t* nrhs, double* b, aocl_int_t* ldb, double* bx, aocl_int_t* ldbx, double* u, aocl_int_t* ldu, double* vt, aocl_int_t* k, double* difl, double* difr, double* z, double* poles, aocl_int_t* givptr, aocl_int_t* givcol, aocl_int_t* ldgcol, aocl_int_t* perm, double* givnum, double* c, double* s, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dlalsd(char* uplo, aocl_int_t* smlsiz, aocl_int_t* n, aocl_int_t* nrhs, double* d, double* e, double* b, aocl_int_t* ldb, double* rcond, aocl_int_t* rank, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dlamrg(aocl_int_t* n1, aocl_int_t* n2, double* a, aocl_int_t* dtrd1, aocl_int_t* dtrd2, aocl_int_t* index);
void LAPACK_EXPORT_dlamswlq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* mb, aocl_int_t* nb, double* a, aocl_int_t* lda, double* t, aocl_int_t* ldt, double* c, aocl_int_t* ldc, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dlamtsqr(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* mb, aocl_int_t* nb, double* a, aocl_int_t* lda, double* t, aocl_int_t* ldt, double* c, aocl_int_t* ldc, double* work, aocl_int_t* lwork, aocl_int_t* info);
aocl_int_t LAPACK_EXPORT_dlaneg(aocl_int_t* n, double* d, double* lld, double* sigma, double* pivmin, aocl_int_t* r__);
void LAPACK_EXPORT_dlanv2(double* a, double* b, double* c, double* d, double* rt1r, double* rt1i, double* rt2r, double* rt2i, double* cs, double* sn);
void LAPACK_EXPORT_dlaorhr_col_getrfnp2(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* d, aocl_int_t* info);
void LAPACK_EXPORT_dlaorhr_col_getrfnp(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* d, aocl_int_t* info);
void LAPACK_EXPORT_dlapll(aocl_int_t* n, double* x, aocl_int_t* incx, double* y, aocl_int_t* incy, double* ssmin);
void LAPACK_EXPORT_dlapmr(logical* forwrd, aocl_int_t* m, aocl_int_t* n, double* x, aocl_int_t* ldx, aocl_int_t* k);
void LAPACK_EXPORT_dlapmt(logical* forwrd, aocl_int_t* m, aocl_int_t* n, double* x, aocl_int_t* ldx, aocl_int_t* k);
void LAPACK_EXPORT_dla_porfsx_extended(aocl_int_t* prec_type, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, logical* colequ, double* c, double* b, aocl_int_t* ldb, double* y, aocl_int_t* ldy, double* berr_out, aocl_int_t* n_norms, double* err_bnds_norm, double* err_bnds_comp, double* res, double* ayb, double* dy, double* y_tail, double* rcond, aocl_int_t* ithresh, double* rthresh, double* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_dlaqgb(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, double* ab, aocl_int_t* ldab, double* r, double* c, double* rowcnd, double* colcnd, double* amax, char* equed);
void LAPACK_EXPORT_dlaqge(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, char* equed);
void LAPACK_EXPORT_dlaqp2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* offset, double* a, aocl_int_t* lda, aocl_int_t* jpvt, double* tau, double* vn1, double* vn2, double* work);
void LAPACK_EXPORT_dlaqps(aocl_int_t* m, aocl_int_t* n, aocl_int_t* offset, aocl_int_t* nb, aocl_int_t* kb, double* a, aocl_int_t* lda, aocl_int_t* jpvt, double* tau, double* vn1, double* vn2, double* auxv, double* f, aocl_int_t* ldf);
void LAPACK_EXPORT_dlaqr0(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, double* h, aocl_int_t* ldh, double* wr, double* wi, aocl_int_t* iloz, aocl_int_t* ihiz, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dlaqr1(aocl_int_t* n, double* h, aocl_int_t* ldh, double* sr1, double* si1, double* sr2, double* si2, double* v);
void LAPACK_EXPORT_dlaqr2(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ktop, aocl_int_t* kbot, aocl_int_t* nw, double* h, aocl_int_t* ldh, aocl_int_t* iloz, aocl_int_t* ihiz, double* z, aocl_int_t* ldz, aocl_int_t* ns, aocl_int_t* nd, double* sr, double* si, double* v, aocl_int_t* ldv, aocl_int_t* nh, double* t, aocl_int_t* ldt, aocl_int_t* nv, double* wv, aocl_int_t* ldwv, double* work, aocl_int_t* lwork);
void LAPACK_EXPORT_dlaqr3(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ktop, aocl_int_t* kbot, aocl_int_t* nw, double* h, aocl_int_t* ldh, aocl_int_t* iloz, aocl_int_t* ihiz, double* z, aocl_int_t* ldz, aocl_int_t* ns, aocl_int_t* nd, double* sr, double* si, double* v, aocl_int_t* ldv, aocl_int_t* nh, double* t, aocl_int_t* ldt, aocl_int_t* nv, double* wv, aocl_int_t* ldwv, double* work, aocl_int_t* lwork);
void LAPACK_EXPORT_dlaqr4(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, double* h, aocl_int_t* ldh, double* wr, double* wi, aocl_int_t* iloz, aocl_int_t* ihiz, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dlaqr5(logical* wantt, logical* wantz, aocl_int_t* kacc22, aocl_int_t* n, aocl_int_t* ktop, aocl_int_t* kbot, aocl_int_t* nshfts, double* sr, double* si, double* h, aocl_int_t* ldh, aocl_int_t* iloz, aocl_int_t* ihiz, double* z, aocl_int_t* ldz, double* v, aocl_int_t* ldv, double* u, aocl_int_t* ldu, aocl_int_t* nv, double* wv, aocl_int_t* ldwv, aocl_int_t* nh, double* wh, aocl_int_t* ldwh);
void LAPACK_EXPORT_dlaqsb(char* uplo, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, double* s, double* scond, double* amax, char* equed);
void LAPACK_EXPORT_dlaqsp(char* uplo, aocl_int_t* n, double* ap, double* s, double* scond, double* amax, char* equed);
void LAPACK_EXPORT_dlaqsy(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* s, double* scond, double* amax, char* equed);
void LAPACK_EXPORT_dlaqtr(logical* ltran, logical* lfloat, aocl_int_t* n, double* t, aocl_int_t* ldt, double* b, double* w, double* scale, double* x, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dlar1v(aocl_int_t* n, aocl_int_t* b1, aocl_int_t* bn, double* lambda, double* d, double* l, double* ld, double* lld, double* pivmin, double* gaptol, double* z, logical* wantnc, aocl_int_t* negcnt, double* ztz, double* mingma, aocl_int_t* r, aocl_int_t* isuppz, double* nrminv, double* resid, double* rqcorr, double* work);
void LAPACK_EXPORT_dlar2v(aocl_int_t* n, double* x, double* y, double* z, aocl_int_t* incx, double* c, double* s, aocl_int_t* incc);
void LAPACK_EXPORT_dlarfb(char* side, char* trans, char* direct, char* storev, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* v, aocl_int_t* ldv, double* t, aocl_int_t* ldt, double* c, aocl_int_t* ldc, double* work, aocl_int_t* ldwork);
void LAPACK_EXPORT_dlarf(char* side, aocl_int_t* m, aocl_int_t* n, double* v, aocl_int_t* incv, double* tau, double* c, aocl_int_t* ldc, double* work);
void LAPACK_EXPORT_dlarfg(aocl_int_t* n, double* alpha, double* x, aocl_int_t* incx, double* tau);
void LAPACK_EXPORT_dlarfgp(aocl_int_t* n, double* alpha, double* x, aocl_int_t* incx, double* tau);
void LAPACK_EXPORT_dlarft(char* direct, char* storev, aocl_int_t* n, aocl_int_t* k, double* v, aocl_int_t* ldv, double* tau, double* t, aocl_int_t* ldt);
void LAPACK_EXPORT_dlarfx(char* side, aocl_int_t* m, aocl_int_t* n, double* v, double* tau, double* c, aocl_int_t* ldc, double* work);
void LAPACK_EXPORT_dlarfy(char* uplo, aocl_int_t* n, double* v, aocl_int_t* incv, double* tau, double* c, aocl_int_t* ldc, double* work);
void LAPACK_EXPORT_dlargv(aocl_int_t* n, double* x, aocl_int_t* incx, double* y, aocl_int_t* incy, double* c, aocl_int_t* incc);
void LAPACK_EXPORT_dlarnv(aocl_int_t* idist, aocl_int_t* iseed, aocl_int_t* n, double* x);
void LAPACK_EXPORT_dlarra(aocl_int_t* n, double* d, double* e, double* e2, double* spltol, double* tnrm, aocl_int_t* nsplit, aocl_int_t* isplit, aocl_int_t* info);
void LAPACK_EXPORT_dlarrb(aocl_int_t* n, double* d, double* lld, aocl_int_t* ifirst, aocl_int_t* ilast, double* rtol1, double* rtol2, aocl_int_t* offset, double* w, double* wgap, double* werr, double* work, aocl_int_t* iwork, double* pivmin, double* spdiam, aocl_int_t* twist, aocl_int_t* info);
void LAPACK_EXPORT_dlarrc(char* jobt, aocl_int_t* n, double* vl, double* vu, double* d, double* e, double* pivmin, aocl_int_t* eigcnt, aocl_int_t* lcnt, aocl_int_t* rcnt, aocl_int_t* info);
void LAPACK_EXPORT_dlarrd(char* range, char* order, aocl_int_t* n, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* gers, double* reltol, double* d, double* e, double* e2, double* pivmin, aocl_int_t* nsplit, aocl_int_t* isplit, aocl_int_t* m, double* w, double* werr, double* wl, double* wu, aocl_int_t* iblock, aocl_int_t* indexw, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dlarre(char* range, aocl_int_t* n, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* d, double* e, double* e2, double* rtol1, double* rtol2, double* spltol, aocl_int_t* nsplit, aocl_int_t* isplit, aocl_int_t* m, double* w, double* werr, double* wgap, aocl_int_t* iblock, aocl_int_t* indexw, double* gers, double* pivmin, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dlarrf(aocl_int_t* n, double* d, double* l, double* ld, aocl_int_t* clstrt, aocl_int_t* clend, double* w, double* wgap, double* werr, double* spdiam, double* clgapl, double* clgapr, double* pivmin, double* sigma, double* dplus, double* lplus, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dlarrj(aocl_int_t* n, double* d, double* e2, aocl_int_t* ifirst, aocl_int_t* ilast, double* rtol, aocl_int_t* offset, double* w, double* werr, double* work, aocl_int_t* iwork, double* pivmin, double* spdiam, aocl_int_t* info);
void LAPACK_EXPORT_dlarrk(aocl_int_t* n, aocl_int_t* iw, double* gl, double* gu, double* d, double* e2, double* pivmin, double* reltol, double* w, double* werr, aocl_int_t* info);
void LAPACK_EXPORT_dlarrr(aocl_int_t* n, double* d, double* e, aocl_int_t* info);
void LAPACK_EXPORT_dlarrv(aocl_int_t* n, double* vl, double* vu, double* d, double* l, double* pivmin, aocl_int_t* isplit, aocl_int_t* m, aocl_int_t* dol, aocl_int_t* dou, double* minrgp, double* rtol1, double* rtol2, double* w, double* werr, double* wgap, aocl_int_t* iblock, aocl_int_t* indexw, double* gers, double* z, aocl_int_t* ldz, aocl_int_t* isuppz, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dlarscl2(aocl_int_t* m, aocl_int_t* n, double* d, double* x, aocl_int_t* ldx);
void LAPACK_EXPORT_dlartg(double* f, double* g, double* cs, double* sn, double* r__);
void LAPACK_EXPORT_dlartgp(double* f, double* g, double* cs, double* sn, double* r__);
void LAPACK_EXPORT_dlartgs(double* x, double* y, double* sigma, double* cs, double* sn);
void LAPACK_EXPORT_dlartv(aocl_int_t* n, double* x, aocl_int_t* incx, double* y, aocl_int_t* incy, double* c, double* s, aocl_int_t* incc);
void LAPACK_EXPORT_dlaruv(aocl_int_t* iseed, aocl_int_t* n, double* x);
void LAPACK_EXPORT_dlarzb(char* side, char* trans, char* direct, char* storev, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, double* v, aocl_int_t* ldv, double* t, aocl_int_t* ldt, double* c, aocl_int_t* ldc, double* work, aocl_int_t* ldwork);
void LAPACK_EXPORT_dlarz(char* side, aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, double* v, aocl_int_t* incv, double* tau, double* c, aocl_int_t* ldc, double* work);
void LAPACK_EXPORT_dlarzt(char* direct, char* storev, aocl_int_t* n, aocl_int_t* k, double* v, aocl_int_t* ldv, double* tau, double* t, aocl_int_t* ldt);
void LAPACK_EXPORT_dlas2(double* f, double* g, double* h, double* ssmin, double* ssmax);
void LAPACK_EXPORT_dlascl2(aocl_int_t* m, aocl_int_t* n, double* d, double* x, aocl_int_t* ldx);
void LAPACK_EXPORT_dlascl(char* type, aocl_int_t* kl, aocl_int_t* ku, double* cfrom, double* cto, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_dlasd0(aocl_int_t* n, aocl_int_t* sqre, double* d, double* e, double* u, aocl_int_t* ldu, double* vt, aocl_int_t* ldvt, aocl_int_t* smlsiz, aocl_int_t* iwork, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dlasd1(aocl_int_t* nl, aocl_int_t* nr, aocl_int_t* sqre, double* d, double* alpha, double* beta, double* u, aocl_int_t* ldu, double* vt, aocl_int_t* ldvt, aocl_int_t* idxq, aocl_int_t* iwork, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dlasd2(aocl_int_t* nl, aocl_int_t* nr, aocl_int_t* sqre, aocl_int_t* k, double* d, double* z, double* alpha, double* beta, double* u, aocl_int_t* ldu, double* vt, aocl_int_t* ldvt, double* dsigma, double* u2, aocl_int_t* ldu2, double* vt2, aocl_int_t* ldvt2, aocl_int_t* idxp, aocl_int_t* idx, aocl_int_t* idxc, aocl_int_t* idxq, aocl_int_t* coltyp, aocl_int_t* info);
void LAPACK_EXPORT_dlasd3(aocl_int_t* nl, aocl_int_t* nr, aocl_int_t* sqre, aocl_int_t* k, double* d, double* q, aocl_int_t* ldq, double* dsigma, double* u, aocl_int_t* ldu, double* u2, aocl_int_t* ldu2, double* vt, aocl_int_t* ldvt, double* vt2, aocl_int_t* ldvt2, aocl_int_t* idxc, aocl_int_t* ctot, double* z, aocl_int_t* info);
void LAPACK_EXPORT_dlasd4(aocl_int_t* n, aocl_int_t* i, double* d, double* z, double* delta, double* rho, double* sigma, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dlasd5(aocl_int_t* i, double* d, double* z, double* delta, double* rho, double* dsigma, double* work);
void LAPACK_EXPORT_dlasd6(aocl_int_t* icompq, aocl_int_t* nl, aocl_int_t* nr, aocl_int_t* sqre, double* d, double* vf, double* vl, double* alpha, double* beta, aocl_int_t* idxq, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, aocl_int_t* ldgcol, double* givnum, aocl_int_t* ldgnum, double* poles, double* difl, double* difr, double* z, aocl_int_t* k, double* c, double* s, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dlasd7(aocl_int_t* icompq, aocl_int_t* nl, aocl_int_t* nr, aocl_int_t* sqre, aocl_int_t* k, double* d, double* z, double* zw, double* vf, double* vfw, double* vl, double* vlw, double* alpha, double* beta, double* dsigma, aocl_int_t* idx, aocl_int_t* idxp, aocl_int_t* idxq, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, aocl_int_t* ldgcol, double* givnum, aocl_int_t* ldgnum, double* c, double* s, aocl_int_t* info);
void LAPACK_EXPORT_dlasd8(aocl_int_t* icompq, aocl_int_t* k, double* d, double* z, double* vf, double* vl, double* difl, double* difr, aocl_int_t* lddifr, double* dsigma, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dlasda(aocl_int_t* icompq, aocl_int_t* smlsiz, aocl_int_t* n, aocl_int_t* sqre, double* d, double* e, double* u, aocl_int_t* ldu, double* vt, aocl_int_t* k, double* difl, double* difr, double* z, double* poles, aocl_int_t* givptr, aocl_int_t* givcol, aocl_int_t* ldgcol, aocl_int_t* perm, double* givnum, double* c, double* s, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dlasdq(char* uplo, aocl_int_t* sqre, aocl_int_t* n, aocl_int_t* ncvt, aocl_int_t* nru, aocl_int_t* ncc, double* d, double* e, double* vt, aocl_int_t* ldvt, double* u, aocl_int_t* ldu, double* c, aocl_int_t* ldc, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dlasdt(aocl_int_t* n, aocl_int_t* lvl, aocl_int_t* nd, aocl_int_t* inode, aocl_int_t* ndiml, aocl_int_t* ndimr, aocl_int_t* msub);
void LAPACK_EXPORT_dlaset(char* uplo, aocl_int_t* m, aocl_int_t* n, double* alpha, double* beta, double* a, aocl_int_t* lda);
void LAPACK_EXPORT_dlasq1(aocl_int_t* n, double* d, double* e, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dlasq2(aocl_int_t* n, double* z, aocl_int_t* info);
void LAPACK_EXPORT_dlasq3(aocl_int_t* i0, aocl_int_t* n0, double* z, aocl_int_t* pp, double* dmin, double* sigma, double* desig, double* qmax, aocl_int_t* nfail, aocl_int_t* iter, aocl_int_t* ndiv, logical* ieee, aocl_int_t* ttype, double* dmin1, double* dmin2, double* dn, double* dn1, double* dn2, double* g, double* tau);
void LAPACK_EXPORT_dlasq4(aocl_int_t* i0, aocl_int_t* n0, double* z, aocl_int_t* pp, aocl_int_t* n0in, double* dmin, double* dmin1, double* dmin2, double* dn, double* dn1, double* dn2, double* tau, aocl_int_t* ttype, double* g);
void LAPACK_EXPORT_dlasq5(aocl_int_t* i0, aocl_int_t* n0, double* z, aocl_int_t* pp, double* tau, double* sigma, double* dmin, double* dmin1, double* dmin2, double* dn, double* dnm1, double* dnm2, logical* ieee, double* eps);
void LAPACK_EXPORT_dlasq6(aocl_int_t* i0, aocl_int_t* n0, double* z, aocl_int_t* pp, double* dmin, double* dmin1, double* dmin2, double* dn, double* dnm1, double* dnm2);
void LAPACK_EXPORT_dlasr(char* side, char* pivot, char* direct, aocl_int_t* m, aocl_int_t* n, double* c, double* s, double* a, aocl_int_t* lda);
void LAPACK_EXPORT_dlasrt(char* id, aocl_int_t* n, double* d, aocl_int_t* info);
void LAPACK_EXPORT_dlassq(aocl_int_t* n, double* x, aocl_int_t* incx, double* scale, double* sumsq);
void LAPACK_EXPORT_dlasv2(double* f, double* g, double* h, double* ssmin, double* ssmax, double* snr, double* csr, double* snl, double* csl);
void LAPACK_EXPORT_dlaswlq(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, aocl_int_t* nb, double* a, aocl_int_t* lda, double* t, aocl_int_t* ldt, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dlaswp(aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* k1, aocl_int_t* k2, aocl_int_t* ipiv, aocl_int_t* incx);
void LAPACK_EXPORT_dlasy2(logical* ltranl, logical* ltranr, aocl_int_t* isgn, aocl_int_t* n1, aocl_int_t* n2, double* tl, aocl_int_t* ldtl, double* tr, aocl_int_t* ldtr, double* b, aocl_int_t* ldb, double* scale, double* x, aocl_int_t* ldx, double* xnorm, aocl_int_t* info);
void LAPACK_EXPORT_dla_syamv(aocl_int_t* uplo, aocl_int_t* n, double* alpha, double* a, aocl_int_t* lda, double* x, aocl_int_t* incx, double* beta, double* y, aocl_int_t* incy);
void LAPACK_EXPORT_dlasyf_aa(char* uplo, aocl_int_t* j1, aocl_int_t* m, aocl_int_t* nb, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* h, aocl_int_t* ldh, double* work);
void LAPACK_EXPORT_dlasyf(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_dlasyf_rk(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, double* a, aocl_int_t* lda, double* e, aocl_int_t* ipiv, double* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_dlasyf_rook(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_dla_syrfsx_extended(aocl_int_t* prec_type, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, aocl_int_t* ipiv, logical* colequ, double* c, double* b, aocl_int_t* ldb, double* y, aocl_int_t* ldy, double* berr_out, aocl_int_t* n_norms, double* err_bnds_norm, double* err_bnds_comp, double* res, double* ayb, double* dy, double* y_tail, double* rcond, aocl_int_t* ithresh, double* rthresh, double* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_dlat2s(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, float* sa, aocl_int_t* ldsa, aocl_int_t* info);
void LAPACK_EXPORT_dlatbs(char* uplo, char* trans, char* diag, char* normin, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, double* x, double* scale, double* cnorm, aocl_int_t* info);
void LAPACK_EXPORT_dlatdf(aocl_int_t* ijob, aocl_int_t* n, double* z, aocl_int_t* ldz, double* rhs, double* rdsum, double* rdscal, aocl_int_t* ipiv, aocl_int_t* jpiv);
void LAPACK_EXPORT_dlatps(char* uplo, char* trans, char* diag, char* normin, aocl_int_t* n, double* ap, double* x, double* scale, double* cnorm, aocl_int_t* info);
void LAPACK_EXPORT_dlatrd(char* uplo, aocl_int_t* n, aocl_int_t* nb, double* a, aocl_int_t* lda, double* e, double* tau, double* w, aocl_int_t* ldw);
void LAPACK_EXPORT_dlatrs(char* uplo, char* trans, char* diag, char* normin, aocl_int_t* n, double* a, aocl_int_t* lda, double* x, double* scale, double* cnorm, aocl_int_t* info);
void LAPACK_EXPORT_dlatrz(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, double* a, aocl_int_t* lda, double* tau, double* work);
void LAPACK_EXPORT_dlatsqr(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, aocl_int_t* nb, double* a, aocl_int_t* lda, double* t, aocl_int_t* ldt, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dlatzm(char* side, aocl_int_t* m, aocl_int_t* n, double* v, aocl_int_t* incv, double* tau, double* c1, double* c2, aocl_int_t* ldc, double* work);
void LAPACK_EXPORT_dlauu2(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_dlauum(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_dla_wwaddw(aocl_int_t* n, double* x, double* y, double* w);
void LAPACK_EXPORT_dopmtr(char* side, char* uplo, char* trans, aocl_int_t* m, aocl_int_t* n, double* ap, double* tau, double* c, aocl_int_t* ldc, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dorbdb1(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, double* x11, aocl_int_t* ldx11, double* x21, aocl_int_t* ldx21, double* theta, double* phi, double* taup1, double* taup2, double* tauq1, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorbdb2(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, double* x11, aocl_int_t* ldx11, double* x21, aocl_int_t* ldx21, double* theta, double* phi, double* taup1, double* taup2, double* tauq1, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorbdb3(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, double* x11, aocl_int_t* ldx11, double* x21, aocl_int_t* ldx21, double* theta, double* phi, double* taup1, double* taup2, double* tauq1, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorbdb4(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, double* x11, aocl_int_t* ldx11, double* x21, aocl_int_t* ldx21, double* theta, double* phi, double* taup1, double* taup2, double* tauq1, double* phantom, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorbdb5(aocl_int_t* m1, aocl_int_t* m2, aocl_int_t* n, double* x1, aocl_int_t* incx1, double* x2, aocl_int_t* incx2, double* q1, aocl_int_t* ldq1, double* q2, aocl_int_t* ldq2, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorbdb6(aocl_int_t* m1, aocl_int_t* m2, aocl_int_t* n, double* x1, aocl_int_t* incx1, double* x2, aocl_int_t* incx2, double* q1, aocl_int_t* ldq1, double* q2, aocl_int_t* ldq2, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorbdb(char* trans, char* signs, aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, double* x11, aocl_int_t* ldx11, double* x12, aocl_int_t* ldx12, double* x21, aocl_int_t* ldx21, double* x22, aocl_int_t* ldx22, double* theta, double* phi, double* taup1, double* taup2, double* tauq1, double* tauq2, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorg2l(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dorg2r(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dorgl2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dorgbr(char* vect, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorglq(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorgql(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorgqr(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorgr2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dorgrq(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorgtr(char* uplo, aocl_int_t* m, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorgtsqr(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, aocl_int_t* nb, double* a, aocl_int_t* lda, double* t, aocl_int_t* ldt, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorhr_col(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nb, double* a, aocl_int_t* lda, double* t, aocl_int_t* ldt, double* d, aocl_int_t* info);
void LAPACK_EXPORT_dorm22(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* n1, aocl_int_t* n2, double* q, aocl_int_t* ldq, double* c, aocl_int_t* ldc, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorm2l(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* c, aocl_int_t* ldc, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dorm2r(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* c, aocl_int_t* ldc, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dormbr(char* vect, char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* c, aocl_int_t* ldc, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dorml2(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* c, aocl_int_t* ldc, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dormlq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* c, aocl_int_t* ldc, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dormql(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* c, aocl_int_t* ldc, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dormqr(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* c, aocl_int_t* ldc, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dormr2(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* c, aocl_int_t* ldc, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dormr3(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, double* a, aocl_int_t* lda, double* tau, double* c, aocl_int_t* ldc, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dormrq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, double* a, aocl_int_t* lda, double* tau, double* c, aocl_int_t* ldc, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dormrz(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, double* a, aocl_int_t* lda, double* tau, double* c, aocl_int_t* ldc, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dormtr(char* side, char* uplo, char* trans, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* tau, double* c, aocl_int_t* ldc, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dpbcon(char* uplo, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, double* anorm, double* rcond, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dpbequ(char* uplo, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, double* s, double* scond, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_dpbrfs(char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, double* ab, aocl_int_t* ldab, double* afb, aocl_int_t* ldafb, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dpbstf(char* uplo, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, aocl_int_t* info);
void LAPACK_EXPORT_dpbsv(char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, double* ab, aocl_int_t* ldab, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dpbsvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, double* ab, aocl_int_t* ldab, double* afb, aocl_int_t* ldafb, char* equed, double* s, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dpbtf2(char* uplo, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, aocl_int_t* info);
void LAPACK_EXPORT_dpbtrf(char* uplo, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, aocl_int_t* info);
void LAPACK_EXPORT_dpbtrs(char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, double* ab, aocl_int_t* ldab, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dpftrf(char* transr, char* uplo, aocl_int_t* n, double* a, aocl_int_t* info);
void LAPACK_EXPORT_dpftri(char* transr, char* uplo, aocl_int_t* n, double* a, aocl_int_t* info);
void LAPACK_EXPORT_dpftrs(char* transr, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dpocon(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* anorm, double* rcond, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dpoequb(aocl_int_t* n, double* a, aocl_int_t* lda, double* s, double* scond, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_dpoequ(aocl_int_t* n, double* a, aocl_int_t* lda, double* s, double* scond, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_dporfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dporfsx(char* uplo, char* equed, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, double* s, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dposv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dposvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, char* equed, double* s, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dposvxx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, char* equed, double* s, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* rpvgrw, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dpotf2(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_dpotrf2(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_dpotrf(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_dpotri(char* uplo, aocl_int_t* n, double* buff_A, aocl_int_t* ldim_A, aocl_int_t* info);
void LAPACK_EXPORT_dpotrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dppcon(char* uplo, aocl_int_t* n, double* ap, double* anorm, double* rcond, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dppequ(char* uplo, aocl_int_t* n, double* ap, double* s, double* scond, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_dpprfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* ap, double* afp, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dppsv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* ap, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dppsvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* ap, double* afp, char* equed, double* s, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dpptrf(char* uplo, aocl_int_t* n, double* ap, aocl_int_t* info);
void LAPACK_EXPORT_dpptri(char* uplo, aocl_int_t* n, double* ap, aocl_int_t* info);
void LAPACK_EXPORT_dpptrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* ap, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dpstf2(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* piv, aocl_int_t* rank, double* tol, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dpstrf(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* piv, aocl_int_t* rank, double* tol, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dptcon(aocl_int_t* n, double* d, double* e, double* anorm, double* rcond, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dpteqr(char* compz, aocl_int_t* n, double* d, double* e, double* z, aocl_int_t* ldz, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dptrfs(aocl_int_t* n, aocl_int_t* nrhs, double* d, double* e, double* df, double* ef, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* ferr, double* berr, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dptsv(aocl_int_t* n, aocl_int_t* nrhs, double* d, double* e, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dptsvx(char* fact, aocl_int_t* n, aocl_int_t* nrhs, double* d, double* e, double* df, double* ef, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dpttrf(aocl_int_t* n, double* d, double* e, aocl_int_t* info);
void LAPACK_EXPORT_dpttrs(aocl_int_t* n, aocl_int_t* nrhs, double* d, double* e, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dptts2(aocl_int_t* n, aocl_int_t* nrhs, double* d, double* e, double* b, aocl_int_t* ldb);
void LAPACK_EXPORT_drscl(aocl_int_t* n, double* sa, double* sx, aocl_int_t* incx);
void LAPACK_EXPORT_dsb2st_kernels(char* uplo, logical* wantz, aocl_int_t* ttype, aocl_int_t* st, aocl_int_t* ed, aocl_int_t* sweep, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* ib, double* a, aocl_int_t* lda, double* v, double* tau, aocl_int_t* ldvt, double* work);
void LAPACK_EXPORT_dsbev_2stage(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsbev(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dsbevd_2stage(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dsbevd(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dsbevx_2stage(char* jobz, char* range, char* uplo, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, double* q, aocl_int_t* ldq, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_dsbevx(char* jobz, char* range, char* uplo, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, double* q, aocl_int_t* ldq, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_dsbgst(char* vect, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, double* ab, aocl_int_t* ldab, double* bb, aocl_int_t* ldbb, double* x, aocl_int_t* ldx, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dsbgv(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, double* ab, aocl_int_t* ldab, double* bb, aocl_int_t* ldbb, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dsbgvd(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, double* ab, aocl_int_t* ldab, double* bb, aocl_int_t* ldbb, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dsbgvx(char* jobz, char* range, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, double* ab, aocl_int_t* ldab, double* bb, aocl_int_t* ldbb, double* q, aocl_int_t* ldq, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_dsbtrd(char* vect, char* uplo, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, double* d, double* e, double* q, aocl_int_t* ldq, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dsfrk(char* transr, char* uplo, char* trans, aocl_int_t* n, aocl_int_t* k, double* alpha, double* a, aocl_int_t* lda, double* beta, double* c__);
void LAPACK_EXPORT_dsgesv(aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* work, float* swork, aocl_int_t* iter, aocl_int_t* info);
void LAPACK_EXPORT_dspcon(char* uplo, aocl_int_t* n, double* ap, aocl_int_t* ipiv, double* anorm, double* rcond, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dspev(char* jobz, char* uplo, aocl_int_t* n, double* ap, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dspevd(char* jobz, char* uplo, aocl_int_t* n, double* ap, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dspevx(char* jobz, char* range, char* uplo, aocl_int_t* n, double* ap, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_dspgst(aocl_int_t* itype, char* uplo, aocl_int_t* n, double* ap, double* bp, aocl_int_t* info);
void LAPACK_EXPORT_dspgvd(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, double* ap, double* bp, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dspgv(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, double* ap, double* bp, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dspgvx(aocl_int_t* itype, char* jobz, char* range, char* uplo, aocl_int_t* n, double* ap, double* bp, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_dsposv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* work, float* swork, aocl_int_t* iter, aocl_int_t* info);
void LAPACK_EXPORT_dsprfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* ap, double* afp, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dspsv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* ap, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dspsvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* ap, double* afp, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dsptrd(char* uplo, aocl_int_t* n, double* ap, double* d, double* e, double* tau, aocl_int_t* info);
void LAPACK_EXPORT_dsptrf(char* uplo, aocl_int_t* n, double* ap, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_dsptri(char* uplo, aocl_int_t* n, double* ap, aocl_int_t* ipiv, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dsptrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* ap, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dstebz(char* range, char* order, aocl_int_t* n, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, double* d, double* e, aocl_int_t* m, aocl_int_t* nsplit, double* w, aocl_int_t* iblock, aocl_int_t* isplit, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dstedc(char* compz, aocl_int_t* n, double* d, double* e, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dstegr(char* jobz, char* range, aocl_int_t* n, double* d, double* e, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, double* z, aocl_int_t* ldz, aocl_int_t* isuppz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dstein(aocl_int_t* n, double* d, double* e, aocl_int_t* m, double* w, aocl_int_t* iblock, aocl_int_t* isplit, double* z, aocl_int_t* ldz, double* work, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_dstemr(char* jobz, char* range, aocl_int_t* n, double* d, double* e, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, aocl_int_t* m, double* w, double* z, aocl_int_t* ldz, aocl_int_t* nzc, aocl_int_t* isuppz, logical* tryrac, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dsteqr(char* compz, aocl_int_t* n, double* d, double* e, double* z, aocl_int_t* ldz, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dsteqr(char* jobz, aocl_int_t* n, double* d, double* e, double* z, aocl_int_t* ldz, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dsterf(aocl_int_t* n, double* d, double* e, aocl_int_t* info);
void LAPACK_EXPORT_dstev(char* jobz, aocl_int_t* n, double* d, double* e, double* z, aocl_int_t* ldz, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dstevd(char* jobz, aocl_int_t* n, double* d, double* e, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dstevr(char* jobz, char* range, aocl_int_t* n, double* d, double* e, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, double* z, aocl_int_t* ldz, aocl_int_t* isuppz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dstevx(char* jobz, char* range, aocl_int_t* n, double* d, double* e, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_dsycon_3(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* e, aocl_int_t* ipiv, double* anorm, double* rcond, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dsycon(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* anorm, double* rcond, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dsycon_rook(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* anorm, double* rcond, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dsyconv(char* uplo, char* way, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dsyconvf(char* uplo, char* way, aocl_int_t* n, double* a, aocl_int_t* lda, double* e, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_dsyconvf_rook(char* uplo, char* way, aocl_int_t* n, double* a, aocl_int_t* lda, double* e, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_dsyequb(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* s, double* scond, double* amax, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dsyev_2stage(char* jobz, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* w, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsyev(char* jobz, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* w, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsyevd_2stage(char* jobz, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* w, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dsyevd(char* jobz, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* w, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dsyevr_2stage(char* jobz, char* range, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, double* z, aocl_int_t* ldz, aocl_int_t* isuppz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dsyevr(char* jobz, char* range, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, double* z, aocl_int_t* ldz, aocl_int_t* isuppz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dsyevx_2stage(char* jobz, char* range, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_dsyevx(char* jobz, char* range, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_dsygs2(aocl_int_t* itype, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dsygst(aocl_int_t* itype, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dsygv_2stage(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* w, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsygvd(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* w, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dsygv(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* w, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsygvx(aocl_int_t* itype, char* jobz, char* range, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, double* z, aocl_int_t* ldz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_dsyrfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dsyrfsx(char* uplo, char* equed, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, aocl_int_t* ipiv, double* s, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dsysv_aa_2stage(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, double* b, aocl_int_t* ldb, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsysv_aa(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsysv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsysv_rk(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* e, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsysv_rook(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsysvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dsysvxx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* af, aocl_int_t* ldaf, aocl_int_t* ipiv, char* equed, double* s, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* rcond, double* rpvgrw, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dsyswapr(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* i1, aocl_int_t* i2);
void LAPACK_EXPORT_dsytd2(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* d, double* e, double* tau, aocl_int_t* info);
void LAPACK_EXPORT_dsytf2(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_dsytf2_rk(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* e, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_dsytf2_rook(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_dsytrd_2stage(char* vect, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* d, double* e, double* tau, double* hous2, aocl_int_t* lhous2, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsytrd(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* d, double* e, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsytrd_sb2st(char* stage1, char* vect, char* uplo, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, double* d, double* e, double* hous, aocl_int_t* lhous, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsytrd_sy2sb(char* uplo, aocl_int_t* n, aocl_int_t* kd, double* a, aocl_int_t* lda, double* ab, aocl_int_t* ldab, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsytrf_aa_2stage(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsytrf_aa(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsytrf(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsytrf_rk(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* e, aocl_int_t* ipiv, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsytrf_rook(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsytri2(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsytri2x(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* work, aocl_int_t* nb, aocl_int_t* info);
void LAPACK_EXPORT_dsytri_3(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* e, aocl_int_t* ipiv, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsytri_3x(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* e, aocl_int_t* ipiv, double* work, aocl_int_t* nb, aocl_int_t* info);
void LAPACK_EXPORT_dsytri(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dsytri_rook(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dsytrs2(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dsytrs_3(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* e, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dsytrs_aa_2stage(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dsytrs_aa(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dsytrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dsytrs_rook(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, aocl_int_t* ipiv, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dtbcon(char* norm, char* uplo, char* diag, aocl_int_t* n, aocl_int_t* kd, double* ab, aocl_int_t* ldab, double* rcond, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dtbrfs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, double* ab, aocl_int_t* ldab, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dtbtrs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, double* ab, aocl_int_t* ldab, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dtfsm(char* transr, char* side, char* uplo, char* trans, char* diag, aocl_int_t* m, aocl_int_t* n, double* alpha, double* a, double* b, aocl_int_t* ldb);
void LAPACK_EXPORT_dtftri(char* transr, char* uplo, char* diag, aocl_int_t* n, double* a, aocl_int_t* info);
void LAPACK_EXPORT_dtfttp(char* transr, char* uplo, aocl_int_t* n, double* arf, double* ap, aocl_int_t* info);
void LAPACK_EXPORT_dtfttr(char* transr, char* uplo, aocl_int_t* n, double* arf, double* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_dtgevc(char* side, char* howmny, logical* select, aocl_int_t* n, double* s, aocl_int_t* lds, double* p, aocl_int_t* ldp, double* vl, aocl_int_t* ldvl, double* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dtgex2(logical* wantq, logical* wantz, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* q, aocl_int_t* ldq, double* z, aocl_int_t* ldz, aocl_int_t* j1, aocl_int_t* n1, aocl_int_t* n2, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dtgexc(logical* wantq, logical* wantz, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* q, aocl_int_t* ldq, double* z, aocl_int_t* ldz, aocl_int_t* ifst, aocl_int_t* ilst, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dtgsen(aocl_int_t* ijob, logical* wantq, logical* wantz, logical* select, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* alphar, double* alphai, double* beta, double* q, aocl_int_t* ldq, double* z, aocl_int_t* ldz, aocl_int_t* m, double* pl, double* pr, double* dif, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dtgsja(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* tola, double* tolb, double* alpha, double* beta, double* u, aocl_int_t* ldu, double* v, aocl_int_t* ldv, double* q, aocl_int_t* ldq, double* work, aocl_int_t* ncycle, aocl_int_t* info);
void LAPACK_EXPORT_dtgsna(char* job, char* howmny, logical* select, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* vl, aocl_int_t* ldvl, double* vr, aocl_int_t* ldvr, double* s, double* dif, aocl_int_t* mm, aocl_int_t* m, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dtgsy2(char* trans, aocl_int_t* ijob, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* c, aocl_int_t* ldc, double* d, aocl_int_t* ldd, double* e, aocl_int_t* lde, double* f, aocl_int_t* ldf, double* scale, double* rdsum, double* rdscal, aocl_int_t* iwork, aocl_int_t* pq, aocl_int_t* info);
void LAPACK_EXPORT_dtgsyl(char* trans, aocl_int_t* ijob, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* c, aocl_int_t* ldc, double* d, aocl_int_t* ldd, double* e, aocl_int_t* lde, double* f, aocl_int_t* ldf, double* scale, double* dif, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dtpcon(char* norm, char* uplo, char* diag, aocl_int_t* n, double* ap, double* rcond, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dtplqt2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_dtplqt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, aocl_int_t* mb, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* t, aocl_int_t* ldt, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dtpmlqt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, aocl_int_t* mb, double* v, aocl_int_t* ldv, double* t, aocl_int_t* ldt, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dtpmqrt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, aocl_int_t* nb, double* v, aocl_int_t* ldv, double* t, aocl_int_t* ldt, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dtpqrt2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_dtpqrt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, aocl_int_t* nb, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* t, aocl_int_t* ldt, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dtprfb(char* side, char* trans, char* direct, char* storev, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, double* v, aocl_int_t* ldv, double* t, aocl_int_t* ldt, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* work, aocl_int_t* ldwork);
void LAPACK_EXPORT_dtprfs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, double* ap, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dtptri(char* uplo, char* diag, aocl_int_t* n, double* ap, aocl_int_t* info);
void LAPACK_EXPORT_dtptrs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, double* ap, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dtpttf(char* transr, char* uplo, aocl_int_t* n, double* ap, double* arf, aocl_int_t* info);
void LAPACK_EXPORT_dtpttr(char* uplo, aocl_int_t* n, double* ap, double* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_dtrcon(char* norm, char* uplo, char* diag, aocl_int_t* n, double* a, aocl_int_t* lda, double* rcond, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dtrevc3(char* side, char* howmny, logical* select, aocl_int_t* n, double* t, aocl_int_t* ldt, double* vl, aocl_int_t* ldvl, double* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, double* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_dtrevc(char* side, char* howmny, logical* select, aocl_int_t* n, double* t, aocl_int_t* ldt, double* vl, aocl_int_t* ldvl, double* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dtrexc(char* compq, aocl_int_t* n, double* t, aocl_int_t* ldt, double* q, aocl_int_t* ldq, aocl_int_t* ifst, aocl_int_t* ilst, double* work, aocl_int_t* info);
void LAPACK_EXPORT_dtrrfs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* x, aocl_int_t* ldx, double* ferr, double* berr, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dtrsen(char* job, char* compq, logical* select, aocl_int_t* n, double* t, aocl_int_t* ldt, double* q, aocl_int_t* ldq, double* wr, double* wi, aocl_int_t* m, double* s, double* sep, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_dtrsna(char* job, char* howmny, logical* select, aocl_int_t* n, double* t, aocl_int_t* ldt, double* vl, aocl_int_t* ldvl, double* vr, aocl_int_t* ldvr, double* s, double* sep, aocl_int_t* mm, aocl_int_t* m, double* work, aocl_int_t* ldwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_dtrsyl(char* trana, char* tranb, aocl_int_t* isgn, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, double* c, aocl_int_t* ldc, double* scale, aocl_int_t* info);
void LAPACK_EXPORT_dtrti2(char* uplo, char* diag, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_dtrtri(char* uplo, char* diag, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_dtrtrs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, double* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_dtrttf(char* transr, char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* arf, aocl_int_t* info);
void LAPACK_EXPORT_dtrttp(char* uplo, aocl_int_t* n, double* a, aocl_int_t* lda, double* ap, aocl_int_t* info);
void LAPACK_EXPORT_dtzrqf(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* tau, aocl_int_t* info);
void LAPACK_EXPORT_dtzrzf(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, double* tau, double* work, aocl_int_t* lwork, aocl_int_t* info);
aocl_int_t LAPACK_EXPORT_icmax1(aocl_int_t* n, scomplex* cx, aocl_int_t* incx);
aocl_int_t LAPACK_EXPORT_ilaclc(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda);
aocl_int_t LAPACK_EXPORT_ilaclr(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda);
aocl_int_t LAPACK_EXPORT_iladiag(char* diag);
aocl_int_t LAPACK_EXPORT_iladlc(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda);
aocl_int_t LAPACK_EXPORT_iladlr(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda);
aocl_int_t LAPACK_EXPORT_ilaenv2stage(aocl_int_t* ispec, char* name, char* opts, aocl_int_t* n1, aocl_int_t* n2, aocl_int_t* n3, aocl_int_t* n4);
aocl_int_t LAPACK_EXPORT_ilaprec(char* prec);
aocl_int_t LAPACK_EXPORT_ilaslc(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda);
aocl_int_t LAPACK_EXPORT_ilaslr(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda);
aocl_int_t LAPACK_EXPORT_ilatrans(char* trans);
aocl_int_t LAPACK_EXPORT_ilauplo(char* uplo);
void LAPACK_EXPORT_ilaver(aocl_int_t* vers_major, aocl_int_t* vers_minor, aocl_int_t* vers_patch__);
aocl_int_t LAPACK_EXPORT_ilazlc(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda);
aocl_int_t LAPACK_EXPORT_ilazlr(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda);
aocl_int_t LAPACK_EXPORT_iparam2stage(aocl_int_t* ispec, char* name, char* opts, aocl_int_t* ni, aocl_int_t* nbi, aocl_int_t* ibi, aocl_int_t* nxi);
aocl_int_t LAPACK_EXPORT_izmax1(aocl_int_t* n, dcomplex* cx, aocl_int_t* incx);
void LAPACK_EXPORT_sbbcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, float* theta, float* phi, float* u1, aocl_int_t* ldu1, float* u2, aocl_int_t* ldu2, float* v1t, aocl_int_t* ldv1t, float* v2t, aocl_int_t* ldv2t, float* b11d, float* b11e, float* b12d, float* b12e, float* b21d, float* b21e, float* b22d, float* b22e, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sbdsdc(char* uplo, char* compq, aocl_int_t* n, float* d, float* e, float* u, aocl_int_t* ldu, float* vt, aocl_int_t* ldvt, float* q, aocl_int_t* iq, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sbdsqr(char* uplo, aocl_int_t* n, aocl_int_t* ncvt, aocl_int_t* nru, aocl_int_t* ncc, float* d, float* e, float* vt, aocl_int_t* ldvt, float* u, aocl_int_t* ldu, float* c, aocl_int_t* ldc, float* rwork, aocl_int_t* info);
void LAPACK_EXPORT_sbdsvdx(char* uplo, char* jobz, char* range, aocl_int_t* n, float* d, float* e, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, aocl_int_t* ns, float* s, float* z, aocl_int_t* ldz, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_scombssq(float* v1, float* v2);
void LAPACK_EXPORT_sdisna(char* job, aocl_int_t* m, aocl_int_t* n, float* d, float* sep, aocl_int_t* info);
void LAPACK_EXPORT_sgbbrd(char* vect, aocl_int_t* m, aocl_int_t* n, aocl_int_t* ncc, aocl_int_t* kl, aocl_int_t* ku, float* ab, aocl_int_t* ldab, float* d, float* e, float* q, aocl_int_t* ldq, float* pt, aocl_int_t* ldpt, float* c, aocl_int_t* ldc, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgbcon(char* norm, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, float* ab, aocl_int_t* ldab, aocl_int_t* ipiv, float* anorm, float* rcond, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgbequb(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, float* ab, aocl_int_t* ldab, float* r, float* c, float* rowcnd, float* colcnd, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_sgbequ(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, float* ab, aocl_int_t* ldab, float* r, float* c, float* rowcnd, float* colcnd, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_sgbrfs(char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, float* ab, aocl_int_t* ldab, float* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgbrfsx(char* trans, char* equed, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, float* ab, aocl_int_t* ldab, float* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, float* r, float* c, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgbsv(aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, float* ab, aocl_int_t* ldab, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_sgbsvx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, float* ab, aocl_int_t* ldab, float* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, char* equed, float* r, float* c, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgbsvxx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, float* ab, aocl_int_t* ldab, float* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, char* equed, float* r, float* c, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* rpvgrw, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgbtf2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, float* ab, aocl_int_t* ldab, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_sgbtrf(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, float* ab, aocl_int_t* ldab, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_sgbtrs(char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, float* ab, aocl_int_t* ldab, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_sgebak(char* job, char* side, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, float* scale, aocl_int_t* m, float* v, aocl_int_t* ldv, aocl_int_t* info);
void LAPACK_EXPORT_sgebal(char* job, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ilo, aocl_int_t* ihi, float* scale, aocl_int_t* info);
void LAPACK_EXPORT_sgebd2(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* d, float* e, float* tauq, float* taup, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgebrd(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* d, float* e, float* tauq, float* taup, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgecon(char* norm, aocl_int_t* n, float* a, aocl_int_t* lda, float* anorm, float* rcond, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgeequb(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_sgeequ(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_sgees(char* jobvs, char* sort, L_fp select, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* sdim, float* wr, float* wi, float* vs, aocl_int_t* ldvs, float* work, aocl_int_t* lwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_sgeesx(char* jobvs, char* sort, L_fp select, char* sense, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* sdim, float* wr, float* wi, float* vs, aocl_int_t* ldvs, float* rconde, float* rcondv, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_sgeev(char* jobvl, char* jobvr, aocl_int_t* n, float* a, aocl_int_t* lda, float* wr, float* wi, float* vl, aocl_int_t* ldvl, float* vr, aocl_int_t* ldvr, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgeevx(char* balanc, char* jobvl, char* jobvr, char* sense, aocl_int_t* n, float* a, aocl_int_t* lda, float* wr, float* wi, float* vl, aocl_int_t* ldvl, float* vr, aocl_int_t* ldvr, aocl_int_t* ilo, aocl_int_t* ihi, float* scale, float* abnrm, float* rconde, float* rcondv, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgegs(char* jobvsl, char* jobvsr, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* alphar, float* alphai, float* beta, float* vsl, aocl_int_t* ldvsl, float* vsr, aocl_int_t* ldvsr, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgegv(char* jobvl, char* jobvr, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* alphar, float* alphai, float* beta, float* vl, aocl_int_t* ldvl, float* vr, aocl_int_t* ldvr, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgehd2(aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgehrd(aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgejsv(char* joba, char* jobu, char* jobv, char* jobr, char* jobt, char* jobp, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* sva, float* u, aocl_int_t* ldu, float* v, aocl_int_t* ldv, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgelq2(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgelqf(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgelq(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* t, aocl_int_t* tsize, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgelqt3(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_sgelqt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, float* a, aocl_int_t* lda, float* t, aocl_int_t* ldt, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgels(char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgelsd(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* s, float* rcond, aocl_int_t* rank, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgelss(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* s, float* rcond, aocl_int_t* rank, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgelsx(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, aocl_int_t* jpvt, float* rcond, aocl_int_t* rank, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgelsy(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, aocl_int_t* jpvt, float* rcond, aocl_int_t* rank, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgemlq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* t, aocl_int_t* tsize, float* c, aocl_int_t* ldc, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgemlqt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* mb, float* v, aocl_int_t* ldv, float* t, aocl_int_t* ldt, float* c, aocl_int_t* ldc, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgemqr(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* t, aocl_int_t* tsize, float* c, aocl_int_t* ldc, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgemqrt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* nb, float* v, aocl_int_t* ldv, float* t, aocl_int_t* ldt, float* c, aocl_int_t* ldc, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgeql2(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgeqlf(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgeqp3(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* jpvt, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgeqpf(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* jpvt, float* tau, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgeqr2(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgeqr2p(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgeqrf(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgeqrfp(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgeqr(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* t, aocl_int_t* tsize, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgeqrt2(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_sgeqrt3(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_sgeqrt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nb, float* a, aocl_int_t* lda, float* t, aocl_int_t* ldt, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgerfs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgerfsx(char* trans, char* equed, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, aocl_int_t* ipiv, float* r, float* c, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgerq2(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgerqf(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgesc2(aocl_int_t* n, float* a, aocl_int_t* lda, float* rhs, aocl_int_t* ipiv, aocl_int_t* jpiv, float* scale);
void LAPACK_EXPORT_sgesdd(char* jobz, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* s, float* u, aocl_int_t* ldu, float* vt, aocl_int_t* ldvt, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgesvd(char* jobu, char* jobv, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* s, float* u, aocl_int_t* ldu, float* vt, aocl_int_t* ldvt, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgesvdq(char* joba, char* jobp, char* jobr, char* jobu, char* jobv, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* s, float* u, aocl_int_t* ldu, float* v, aocl_int_t* ldv, aocl_int_t* numrank, aocl_int_t* iwork, aocl_int_t* liwork, float* work, aocl_int_t* lwork, float* rwork, aocl_int_t* lrwork, aocl_int_t* info);
void LAPACK_EXPORT_sgesvdx(char* jobu, char* jobvt, char* range, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, aocl_int_t* ns, float* s, float* u, aocl_int_t* ldu, float* vt, aocl_int_t* ldvt, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgesv(aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_sgesvj(char* joba, char* jobu, char* jobv, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* sva, aocl_int_t* mv, float* v, aocl_int_t* ldv, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgesvx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, aocl_int_t* ipiv, char* equed, float* r, float* c, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgesvxx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, aocl_int_t* ipiv, char* equed, float* r, float* c, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* rpvgrw, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgetc2(aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* jpiv, aocl_int_t* info);
void LAPACK_EXPORT_sgetf2(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_sgetrf2(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_sgetrf(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_sgetri(aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgetrs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_sgetsls(char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sggbak(char* job, char* side, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, float* lscale, float* rscale, aocl_int_t* m, float* v, aocl_int_t* ldv, aocl_int_t* info);
void LAPACK_EXPORT_sggbal(char* job, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, aocl_int_t* ilo, aocl_int_t* ihi, float* lscale, float* rscale, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgges3(char* jobvsl, char* jobvsr, char* sort, L_fp selctg, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, aocl_int_t* sdim, float* alphar, float* alphai, float* beta, float* vsl, aocl_int_t* ldvsl, float* vsr, aocl_int_t* ldvsr, float* work, aocl_int_t* lwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_sgges(char* jobvsl, char* jobvsr, char* sort, L_fp selctg, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, aocl_int_t* sdim, float* alphar, float* alphai, float* beta, float* vsl, aocl_int_t* ldvsl, float* vsr, aocl_int_t* ldvsr, float* work, aocl_int_t* lwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_sggesx(char* jobvsl, char* jobvsr, char* sort, L_fp selctg, char* sense, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, aocl_int_t* sdim, float* alphar, float* alphai, float* beta, float* vsl, aocl_int_t* ldvsl, float* vsr, aocl_int_t* ldvsr, float* rconde, float* rcondv, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_sggev3(char* jobvl, char* jobvr, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* alphar, float* alphai, float* beta, float* vl, aocl_int_t* ldvl, float* vr, aocl_int_t* ldvr, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sggev(char* jobvl, char* jobvr, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* alphar, float* alphai, float* beta, float* vl, aocl_int_t* ldvl, float* vr, aocl_int_t* ldvr, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sggevx(char* balanc, char* jobvl, char* jobvr, char* sense, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* alphar, float* alphai, float* beta, float* vl, aocl_int_t* ldvl, float* vr, aocl_int_t* ldvr, aocl_int_t* ilo, aocl_int_t* ihi, float* lscale, float* rscale, float* abnrm, float* bbnrm, float* rconde, float* rcondv, float* work, aocl_int_t* lwork, aocl_int_t* iwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_sggglm(aocl_int_t* n, aocl_int_t* m, aocl_int_t* p, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* d, float* x, float* y, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgghd3(char* compq, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* q, aocl_int_t* ldq, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgghrd(char* compq, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* q, aocl_int_t* ldq, float* z, aocl_int_t* ldz, aocl_int_t* info);
void LAPACK_EXPORT_sgglse(aocl_int_t* m, aocl_int_t* n, aocl_int_t* p, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* c, float* d, float* x, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sggqrf(aocl_int_t* n, aocl_int_t* m, aocl_int_t* p, float* a, aocl_int_t* lda, float* taua, float* b, aocl_int_t* ldb, float* taub, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sggrqf(aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, float* a, aocl_int_t* lda, float* taua, float* b, aocl_int_t* ldb, float* taub, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sggsvd3(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* n, aocl_int_t* p, aocl_int_t* k, aocl_int_t* l, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* alpha, float* beta, float* u, aocl_int_t* ldu, float* v, aocl_int_t* ldv, float* q, aocl_int_t* ldq, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sggsvd(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* n, aocl_int_t* p, aocl_int_t* k, aocl_int_t* l, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* alpha, float* beta, float* u, aocl_int_t* ldu, float* v, aocl_int_t* ldv, float* q, aocl_int_t* ldq, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sggsvp3(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* tola, float* tolb, aocl_int_t* k, aocl_int_t* l, float* u, aocl_int_t* ldu, float* v, aocl_int_t* ldv, float* q, aocl_int_t* ldq, aocl_int_t* iwork, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sggsvp(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* tola, float* tolb, aocl_int_t* k, aocl_int_t* l, float* u, aocl_int_t* ldu, float* v, aocl_int_t* ldv, float* q, aocl_int_t* ldq, aocl_int_t* iwork, float* tau, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sgsvj0(char* jobv, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* d, float* sva, aocl_int_t* mv, float* v, aocl_int_t* ldv, float* eps, float* sfmin, float* tol, aocl_int_t* nsweep, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgsvj1(char* jobv, aocl_int_t* m, aocl_int_t* n, aocl_int_t* n1, float* a, aocl_int_t* lda, float* d, float* sva, aocl_int_t* mv, float* v, aocl_int_t* ldv, float* eps, float* sfmin, float* tol, aocl_int_t* nsweep, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sgtcon(char* norm, aocl_int_t* n, float* dl, float* d, float* du, float* du2, aocl_int_t* ipiv, float* anorm, float* rcond, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgtrfs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, float* dl, float* d, float* du, float* dlf, float* df, float* duf, float* du2, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgtsv(aocl_int_t* n, aocl_int_t* nrhs, float* dl, float* d, float* du, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_sgtsvx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* nrhs, float* dl, float* d, float* du, float* dlf, float* df, float* duf, float* du2, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sgttrf(aocl_int_t* n, float* dl, float* d, float* du, float* du2, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_sgttrs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, float* dl, float* d, float* du, float* du2, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_sgtts2(aocl_int_t* itrans, aocl_int_t* n, aocl_int_t* nrhs, float* dl, float* d, float* du, float* du2, aocl_int_t* ipiv, float* b, aocl_int_t* ldb);
void LAPACK_EXPORT_shgeqz(char* job, char* compq, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, float* h, aocl_int_t* ldh, float* t, aocl_int_t* ldt, float* alphar, float* alphai, float* beta, float* q, aocl_int_t* ldq, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_shsein(char* side, char* eigsrc, char* initv, logical* select, aocl_int_t* n, float* h, aocl_int_t* ldh, float* wr, float* wi, float* vl, aocl_int_t* ldvl, float* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, float* work, aocl_int_t* ifaill, aocl_int_t* ifailr, aocl_int_t* info);
void LAPACK_EXPORT_shseqr(char* job, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, float* h, aocl_int_t* ldh, float* wr, float* wi, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_slabad(float* small_, float* large);
void LAPACK_EXPORT_slabrd(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nb, float* a, aocl_int_t* lda, float* d, float* e, float* tauq, float* taup, float* x, aocl_int_t* ldx, float* y, aocl_int_t* ldy);
void LAPACK_EXPORT_slacn2(aocl_int_t* n, float* v, float* x, aocl_int_t* isgn, float* est, aocl_int_t* kase, aocl_int_t* isave);
void LAPACK_EXPORT_slacon(aocl_int_t* n, float* v, float* x, aocl_int_t* isgn, float* est, aocl_int_t* kase);
void LAPACK_EXPORT_slacpy(char* uplo, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb);
void LAPACK_EXPORT_sladiv(float* a, float* b, float* c, float* d, float* p, float* q);
void LAPACK_EXPORT_slae2(float* a, float* b, float* c, float* rt1, float* rt2);
void LAPACK_EXPORT_slaebz(aocl_int_t* ijob, aocl_int_t* nitmax, aocl_int_t* n, aocl_int_t* mmax, aocl_int_t* minp, aocl_int_t* nbmin, float* abstol, float* reltol, float* pivmin, float* d, float* e, float* e2, aocl_int_t* nval, float* ab, float* c, aocl_int_t* mout, aocl_int_t* nab, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_slaed0(aocl_int_t* icompq, aocl_int_t* qsiz, aocl_int_t* n, float* d, float* e, float* q, aocl_int_t* ldq, float* qstore, aocl_int_t* ldqs, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_slaed1(aocl_int_t* n, float* d, float* q, aocl_int_t* ldq, aocl_int_t* indxq, float* rho, aocl_int_t* cutpnt, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_slaed2(aocl_int_t* k, aocl_int_t* n, aocl_int_t* n1, float* d, float* q, aocl_int_t* ldq, aocl_int_t* indxq, float* rho, float* z, float* dlamda, float* w, float* q2, aocl_int_t* indx, aocl_int_t* indxc, aocl_int_t* indxp, aocl_int_t* coltyp, aocl_int_t* info);
void LAPACK_EXPORT_slaed3(aocl_int_t* k, aocl_int_t* n, aocl_int_t* n1, float* d, float* q, aocl_int_t* ldq, float* rho, float* dlamda, float* q2, aocl_int_t* indx, aocl_int_t* ctot, float* w, float* s, aocl_int_t* info);
void LAPACK_EXPORT_slaed4(aocl_int_t* n, aocl_int_t* i, float* d, float* z, float* delta, float* rho, float* dlam, aocl_int_t* info);
void LAPACK_EXPORT_slaed5(aocl_int_t* i, float* d, float* z, float* delta, float* rho, float* dlam);
void LAPACK_EXPORT_slaed6(aocl_int_t* kniter, logical* orgati, float* rho, float* d, float* z, float* finit, float* tau, aocl_int_t* info);
void LAPACK_EXPORT_slaed7(aocl_int_t* icompq, aocl_int_t* n, aocl_int_t* qsiz, aocl_int_t* tlvls, aocl_int_t* curlvl, aocl_int_t* curpbm, float* d, float* q, aocl_int_t* ldq, aocl_int_t* indxq, float* rho, aocl_int_t* cutpnt, float* qstore, aocl_int_t* qptr, aocl_int_t* prmptr, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, float* givnum, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_slaed8(aocl_int_t* icompq, aocl_int_t* k, aocl_int_t* n, aocl_int_t* qsiz, float* d, float* q, aocl_int_t* ldq, aocl_int_t* indxq, float* rho, aocl_int_t* cutpnt, float* z, float* dlamda, float* q2, aocl_int_t* ldq2, float* w, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, float* givnum, aocl_int_t* indxp, aocl_int_t* indx, aocl_int_t* info);
void LAPACK_EXPORT_slaed9(aocl_int_t* k, aocl_int_t* kstart, aocl_int_t* kstop, aocl_int_t* n, float* d, float* q, aocl_int_t* ldq, float* rho, float* dlamda, float* w, float* s, aocl_int_t* lds, aocl_int_t* info);
void LAPACK_EXPORT_slaeda(aocl_int_t* n, aocl_int_t* tlvls, aocl_int_t* curlvl, aocl_int_t* curpbm, aocl_int_t* prmptr, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, float* givnum, float* q, aocl_int_t* qptr, float* z, float* ztemp, aocl_int_t* info);
void LAPACK_EXPORT_slaein(logical* rightv, logical* noinit, aocl_int_t* n, float* h, aocl_int_t* ldh, float* wr, float* wi, float* vr, float* vi, float* b, aocl_int_t* ldb, float* work, float* eps3, float* smlnum, float* bignum, aocl_int_t* info);
void LAPACK_EXPORT_slaev2(float* a, float* b, float* c, float* rt1, float* rt2, float* cs1, float* sn1);
void LAPACK_EXPORT_slaexc(logical* wantq, aocl_int_t* n, float* t, aocl_int_t* ldt, float* q, aocl_int_t* ldq, aocl_int_t* j1, aocl_int_t* n1, aocl_int_t* n2, float* work, aocl_int_t* info);
void LAPACK_EXPORT_slag2d(aocl_int_t* m, aocl_int_t* n, float* sa, aocl_int_t* ldsa, double* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_slag2(float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* safmin, float* scale1, float* scale2, float* wr1, float* wr2, float* wi);
void LAPACK_EXPORT_sla_gbamv(aocl_int_t* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, float* alpha, float* ab, aocl_int_t* ldab, float* x, aocl_int_t* incx, float* beta, float* y, aocl_int_t* incy);
void LAPACK_EXPORT_sla_gbrfsx_extended(aocl_int_t* prec_type, aocl_int_t* trans_type, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, float* ab, aocl_int_t* ldab, float* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, logical* colequ, float* c, float* b, aocl_int_t* ldb, float* y, aocl_int_t* ldy, float* berr_out, aocl_int_t* n_norms, float* err_bnds_norm, float* err_bnds_comp, float* res, float* ayb, float* dy, float* y_tail, float* rcond, aocl_int_t* ithresh, float* rthresh, float* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_sla_geamv(aocl_int_t* trans, aocl_int_t* m, aocl_int_t* n, float* alpha, float* a, aocl_int_t* lda, float* x, aocl_int_t* incx, float* beta, float* y, aocl_int_t* incy);
void LAPACK_EXPORT_sla_gerfsx_extended(aocl_int_t* prec_type, aocl_int_t* trans_type, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, aocl_int_t* ipiv, logical* colequ, float* c, float* b, aocl_int_t* ldb, float* y, aocl_int_t* ldy, float* berr_out, aocl_int_t* n_norms, float* errs_n, float* errs_c, float* res, float* ayb, float* dy, float* y_tail, float* rcond, aocl_int_t* ithresh, float* rthresh, float* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_slags2(logical* upper, float* a1, float* a2, float* a3, float* b1, float* b2, float* b3, float* csu, float* snu, float* csv, float* snv, float* csq, float* snq);
void LAPACK_EXPORT_slagtf(aocl_int_t* n, float* a, float* lambda, float* b, float* c, float* tol, float* d, aocl_int_t* in, aocl_int_t* info);
void LAPACK_EXPORT_slagtm(char* trans, aocl_int_t* n, aocl_int_t* nrhs, float* alpha, float* dl, float* d, float* du, float* x, aocl_int_t* ldx, float* beta, float* b, aocl_int_t* ldb);
void LAPACK_EXPORT_slagts(aocl_int_t* job, aocl_int_t* n, float* a, float* b, float* c, float* d, aocl_int_t* in, float* y, float* tol, aocl_int_t* info);
void LAPACK_EXPORT_slagv2(float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* alphar, float* alphai, float* beta, float* csl, float* snl, float* csr, float* snr);
void LAPACK_EXPORT_slahqr(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, float* h, aocl_int_t* ldh, float* wr, float* wi, aocl_int_t* iloz, aocl_int_t* ihiz, float* z, aocl_int_t* ldz, aocl_int_t* info);
void LAPACK_EXPORT_slahr2(aocl_int_t* n, aocl_int_t* k, aocl_int_t* nb, float* a, aocl_int_t* lda, float* tau, float* t, aocl_int_t* ldt, float* y, aocl_int_t* ldy);
void LAPACK_EXPORT_slahrd(aocl_int_t* n, aocl_int_t* k, aocl_int_t* nb, float* a, aocl_int_t* lda, float* tau, float* t, aocl_int_t* ldt, float* y, aocl_int_t* ldy);
void LAPACK_EXPORT_slaic1(aocl_int_t* job, aocl_int_t* j, float* x, float* sest, float* w, float* gamma, float* sestpr, float* s, float* c__);
void LAPACK_EXPORT_sla_lin_berr(aocl_int_t* n, aocl_int_t* nz, aocl_int_t* nrhs, float* res, float* ayb, float* berr);
void LAPACK_EXPORT_slaln2(logical* ltrans, aocl_int_t* na, aocl_int_t* nw, float* smin, float* ca, float* a, aocl_int_t* lda, float* d1, float* d2, float* b, aocl_int_t* ldb, float* wr, float* wi, float* x, aocl_int_t* ldx, float* scale, float* xnorm, aocl_int_t* info);
void LAPACK_EXPORT_slals0(aocl_int_t* icompq, aocl_int_t* nl, aocl_int_t* nr, aocl_int_t* sqre, aocl_int_t* nrhs, float* b, aocl_int_t* ldb, float* bx, aocl_int_t* ldbx, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, aocl_int_t* ldgcol, float* givnum, aocl_int_t* ldgnum, float* poles, float* difl, float* difr, float* z, aocl_int_t* k, float* c, float* s, float* work, aocl_int_t* info);
void LAPACK_EXPORT_slalsa(aocl_int_t* icompq, aocl_int_t* smlsiz, aocl_int_t* n, aocl_int_t* nrhs, float* b, aocl_int_t* ldb, float* bx, aocl_int_t* ldbx, float* u, aocl_int_t* ldu, float* vt, aocl_int_t* k, float* difl, float* difr, float* z, float* poles, aocl_int_t* givptr, aocl_int_t* givcol, aocl_int_t* ldgcol, aocl_int_t* perm, float* givnum, float* c, float* s, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_slalsd(char* uplo, aocl_int_t* smlsiz, aocl_int_t* n, aocl_int_t* nrhs, float* d, float* e, float* b, aocl_int_t* ldb, float* rcond, aocl_int_t* rank, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_slamrg(aocl_int_t* n1, aocl_int_t* n2, float* a, aocl_int_t* strd1, aocl_int_t* strd2, aocl_int_t* index);
void LAPACK_EXPORT_slamswlq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* mb, aocl_int_t* nb, float* a, aocl_int_t* lda, float* t, aocl_int_t* ldt, float* c, aocl_int_t* ldc, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_slamtsqr(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* mb, aocl_int_t* nb, float* a, aocl_int_t* lda, float* t, aocl_int_t* ldt, float* c, aocl_int_t* ldc, float* work, aocl_int_t* lwork, aocl_int_t* info);
aocl_int_t LAPACK_EXPORT_slaneg(aocl_int_t* n, float* d, float* lld, float* sigma, float* pivmin, aocl_int_t* r__);
void LAPACK_EXPORT_slanv2(float* a, float* b, float* c, float* d, float* rt1r, float* rt1i, float* rt2r, float* rt2i, float* cs, float* sn);
void LAPACK_EXPORT_slaorhr_col_getrfnp2(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* d, aocl_int_t* info);
void LAPACK_EXPORT_slaorhr_col_getrfnp(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* d, aocl_int_t* info);
void LAPACK_EXPORT_slapll(aocl_int_t* n, float* x, aocl_int_t* incx, float* y, aocl_int_t* incy, float* ssmin);
void LAPACK_EXPORT_slapmr(logical* forwrd, aocl_int_t* m, aocl_int_t* n, float* x, aocl_int_t* ldx, aocl_int_t* k);
void LAPACK_EXPORT_slapmt(logical* forwrd, aocl_int_t* m, aocl_int_t* n, float* x, aocl_int_t* ldx, aocl_int_t* k);
void LAPACK_EXPORT_sla_porfsx_extended(aocl_int_t* prec_type, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, logical* colequ, float* c, float* b, aocl_int_t* ldb, float* y, aocl_int_t* ldy, float* berr_out, aocl_int_t* n_norms, float* err_bnds_norm, float* err_bnds_comp, float* res, float* ayb, float* dy, float* y_tail, float* rcond, aocl_int_t* ithresh, float* rthresh, float* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_slaqgb(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, float* ab, aocl_int_t* ldab, float* r, float* c, float* rowcnd, float* colcnd, float* amax, char* equed);
void LAPACK_EXPORT_slaqge(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, char* equed);
void LAPACK_EXPORT_slaqp2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* offset, float* a, aocl_int_t* lda, aocl_int_t* jpvt, float* tau, float* vn1, float* vn2, float* work);
void LAPACK_EXPORT_slaqps(aocl_int_t* m, aocl_int_t* n, aocl_int_t* offset, aocl_int_t* nb, aocl_int_t* kb, float* a, aocl_int_t* lda, aocl_int_t* jpvt, float* tau, float* vn1, float* vn2, float* auxv, float* f, aocl_int_t* ldf);
void LAPACK_EXPORT_slaqr0(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, float* h, aocl_int_t* ldh, float* wr, float* wi, aocl_int_t* iloz, aocl_int_t* ihiz, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_slaqr1(aocl_int_t* n, float* h, aocl_int_t* ldh, float* sr1, float* si1, float* sr2, float* si2, float* v);
void LAPACK_EXPORT_slaqr2(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ktop, aocl_int_t* kbot, aocl_int_t* nw, float* h, aocl_int_t* ldh, aocl_int_t* iloz, aocl_int_t* ihiz, float* z, aocl_int_t* ldz, aocl_int_t* ns, aocl_int_t* nd, float* sr, float* si, float* v, aocl_int_t* ldv, aocl_int_t* nh, float* t, aocl_int_t* ldt, aocl_int_t* nv, float* wv, aocl_int_t* ldwv, float* work, aocl_int_t* lwork);
void LAPACK_EXPORT_slaqr3(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ktop, aocl_int_t* kbot, aocl_int_t* nw, float* h, aocl_int_t* ldh, aocl_int_t* iloz, aocl_int_t* ihiz, float* z, aocl_int_t* ldz, aocl_int_t* ns, aocl_int_t* nd, float* sr, float* si, float* v, aocl_int_t* ldv, aocl_int_t* nh, float* t, aocl_int_t* ldt, aocl_int_t* nv, float* wv, aocl_int_t* ldwv, float* work, aocl_int_t* lwork);
void LAPACK_EXPORT_slaqr4(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, float* h, aocl_int_t* ldh, float* wr, float* wi, aocl_int_t* iloz, aocl_int_t* ihiz, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_slaqr5(logical* wantt, logical* wantz, aocl_int_t* kacc22, aocl_int_t* n, aocl_int_t* ktop, aocl_int_t* kbot, aocl_int_t* nshfts, float* sr, float* si, float* h, aocl_int_t* ldh, aocl_int_t* iloz, aocl_int_t* ihiz, float* z, aocl_int_t* ldz, float* v, aocl_int_t* ldv, float* u, aocl_int_t* ldu, aocl_int_t* nv, float* wv, aocl_int_t* ldwv, aocl_int_t* nh, float* wh, aocl_int_t* ldwh);
void LAPACK_EXPORT_slaqsb(char* uplo, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, float* s, float* scond, float* amax, char* equed);
void LAPACK_EXPORT_slaqsp(char* uplo, aocl_int_t* n, float* ap, float* s, float* scond, float* amax, char* equed);
void LAPACK_EXPORT_slaqsy(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* s, float* scond, float* amax, char* equed);
void LAPACK_EXPORT_slaqtr(logical* ltran, logical* lfloat, aocl_int_t* n, float* t, aocl_int_t* ldt, float* b, float* w, float* scale, float* x, float* work, aocl_int_t* info);
void LAPACK_EXPORT_slar1v(aocl_int_t* n, aocl_int_t* b1, aocl_int_t* bn, float* lambda, float* d, float* l, float* ld, float* lld, float* pivmin, float* gaptol, float* z, logical* wantnc, aocl_int_t* negcnt, float* ztz, float* mingma, aocl_int_t* r, aocl_int_t* isuppz, float* nrminv, float* resid, float* rqcorr, float* work);
void LAPACK_EXPORT_slar2v(aocl_int_t* n, float* x, float* y, float* z, aocl_int_t* incx, float* c, float* s, aocl_int_t* incc);
void LAPACK_EXPORT_slarfb(char* side, char* trans, char* direct, char* storev, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* v, aocl_int_t* ldv, float* t, aocl_int_t* ldt, float* c, aocl_int_t* ldc, float* work, aocl_int_t* ldwork);
void LAPACK_EXPORT_slarf(char* side, aocl_int_t* m, aocl_int_t* n, float* v, aocl_int_t* incv, float* tau, float* c, aocl_int_t* ldc, float* work);
void LAPACK_EXPORT_slarfg(aocl_int_t* n, float* alpha, float* x, aocl_int_t* incx, float* tau);
void LAPACK_EXPORT_slarfgp(aocl_int_t* n, float* alpha, float* x, aocl_int_t* incx, float* tau);
void LAPACK_EXPORT_slarft(char* direct, char* storev, aocl_int_t* n, aocl_int_t* k, float* v, aocl_int_t* ldv, float* tau, float* t, aocl_int_t* ldt);
void LAPACK_EXPORT_slarfx(char* side, aocl_int_t* m, aocl_int_t* n, float* v, float* tau, float* c, aocl_int_t* ldc, float* work);
void LAPACK_EXPORT_slarfy(char* uplo, aocl_int_t* n, float* v, aocl_int_t* incv, float* tau, float* c, aocl_int_t* ldc, float* work);
void LAPACK_EXPORT_slargv(aocl_int_t* n, float* x, aocl_int_t* incx, float* y, aocl_int_t* incy, float* c, aocl_int_t* incc);
void LAPACK_EXPORT_slarnv(aocl_int_t* idist, aocl_int_t* iseed, aocl_int_t* n, float* x);
void LAPACK_EXPORT_slarra(aocl_int_t* n, float* d, float* e, float* e2, float* spltol, float* tnrm, aocl_int_t* nsplit, aocl_int_t* isplit, aocl_int_t* info);
void LAPACK_EXPORT_slarrb(aocl_int_t* n, float* d, float* lld, aocl_int_t* ifirst, aocl_int_t* ilast, float* rtol1, float* rtol2, aocl_int_t* offset, float* w, float* wgap, float* werr, float* work, aocl_int_t* iwork, float* pivmin, float* spdiam, aocl_int_t* twist, aocl_int_t* info);
void LAPACK_EXPORT_slarrc(char* jobt, aocl_int_t* n, float* vl, float* vu, float* d, float* e, float* pivmin, aocl_int_t* eigcnt, aocl_int_t* lcnt, aocl_int_t* rcnt, aocl_int_t* info);
void LAPACK_EXPORT_slarrd(char* range, char* order, aocl_int_t* n, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* gers, float* reltol, float* d, float* e, float* e2, float* pivmin, aocl_int_t* nsplit, aocl_int_t* isplit, aocl_int_t* m, float* w, float* werr, float* wl, float* wu, aocl_int_t* iblock, aocl_int_t* indexw, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_slarre(char* range, aocl_int_t* n, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* d, float* e, float* e2, float* rtol1, float* rtol2, float* spltol, aocl_int_t* nsplit, aocl_int_t* isplit, aocl_int_t* m, float* w, float* werr, float* wgap, aocl_int_t* iblock, aocl_int_t* indexw, float* gers, float* pivmin, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_slarrf(aocl_int_t* n, float* d, float* l, float* ld, aocl_int_t* clstrt, aocl_int_t* clend, float* w, float* wgap, float* werr, float* spdiam, float* clgapl, float* clgapr, float* pivmin, float* sigma, float* dplus, float* lplus, float* work, aocl_int_t* info);
void LAPACK_EXPORT_slarrj(aocl_int_t* n, float* d, float* e2, aocl_int_t* ifirst, aocl_int_t* ilast, float* rtol, aocl_int_t* offset, float* w, float* werr, float* work, aocl_int_t* iwork, float* pivmin, float* spdiam, aocl_int_t* info);
void LAPACK_EXPORT_slarrk(aocl_int_t* n, aocl_int_t* iw, float* gl, float* gu, float* d, float* e2, float* pivmin, float* reltol, float* w, float* werr, aocl_int_t* info);
void LAPACK_EXPORT_slarrr(aocl_int_t* n, float* d, float* e, aocl_int_t* info);
void LAPACK_EXPORT_slarrv(aocl_int_t* n, float* vl, float* vu, float* d, float* l, float* pivmin, aocl_int_t* isplit, aocl_int_t* m, aocl_int_t* dol, aocl_int_t* dou, float* minrgp, float* rtol1, float* rtol2, float* w, float* werr, float* wgap, aocl_int_t* iblock, aocl_int_t* indexw, float* gers, float* z, aocl_int_t* ldz, aocl_int_t* isuppz, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_slarscl2(aocl_int_t* m, aocl_int_t* n, float* d, float* x, aocl_int_t* ldx);
void LAPACK_EXPORT_slartg(float* f, float* g, float* cs, float* sn, float* r__);
void LAPACK_EXPORT_slartgp(float* f, float* g, float* cs, float* sn, float* r__);
void LAPACK_EXPORT_slartgs(float* x, float* y, float* sigma, float* cs, float* sn);
void LAPACK_EXPORT_slartv(aocl_int_t* n, float* x, aocl_int_t* incx, float* y, aocl_int_t* incy, float* c, float* s, aocl_int_t* incc);
void LAPACK_EXPORT_slaruv(aocl_int_t* iseed, aocl_int_t* n, float* x);
void LAPACK_EXPORT_slarzb(char* side, char* trans, char* direct, char* storev, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, float* v, aocl_int_t* ldv, float* t, aocl_int_t* ldt, float* c, aocl_int_t* ldc, float* work, aocl_int_t* ldwork);
void LAPACK_EXPORT_slarz(char* side, aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, float* v, aocl_int_t* incv, float* tau, float* c, aocl_int_t* ldc, float* work);
void LAPACK_EXPORT_slarzt(char* direct, char* storev, aocl_int_t* n, aocl_int_t* k, float* v, aocl_int_t* ldv, float* tau, float* t, aocl_int_t* ldt);
void LAPACK_EXPORT_slas2(float* f, float* g, float* h, float* ssmin, float* ssmax);
void LAPACK_EXPORT_slascl2(aocl_int_t* m, aocl_int_t* n, float* d, float* x, aocl_int_t* ldx);
void LAPACK_EXPORT_slascl(char* type, aocl_int_t* kl, aocl_int_t* ku, float* cfrom, float* cto, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_slasd0(aocl_int_t* n, aocl_int_t* sqre, float* d, float* e, float* u, aocl_int_t* ldu, float* vt, aocl_int_t* ldvt, aocl_int_t* smlsiz, aocl_int_t* iwork, float* work, aocl_int_t* info);
void LAPACK_EXPORT_slasd1(aocl_int_t* nl, aocl_int_t* nr, aocl_int_t* sqre, float* d, float* alpha, float* beta, float* u, aocl_int_t* ldu, float* vt, aocl_int_t* ldvt, aocl_int_t* idxq, aocl_int_t* iwork, float* work, aocl_int_t* info);
void LAPACK_EXPORT_slasd2(aocl_int_t* nl, aocl_int_t* nr, aocl_int_t* sqre, aocl_int_t* k, float* d, float* z, float* alpha, float* beta, float* u, aocl_int_t* ldu, float* vt, aocl_int_t* ldvt, float* dsigma, float* u2, aocl_int_t* ldu2, float* vt2, aocl_int_t* ldvt2, aocl_int_t* idxp, aocl_int_t* idx, aocl_int_t* idxc, aocl_int_t* idxq, aocl_int_t* coltyp, aocl_int_t* info);
void LAPACK_EXPORT_slasd3(aocl_int_t* nl, aocl_int_t* nr, aocl_int_t* sqre, aocl_int_t* k, float* d, float* q, aocl_int_t* ldq, float* dsigma, float* u, aocl_int_t* ldu, float* u2, aocl_int_t* ldu2, float* vt, aocl_int_t* ldvt, float* vt2, aocl_int_t* ldvt2, aocl_int_t* idxc, aocl_int_t* ctot, float* z, aocl_int_t* info);
void LAPACK_EXPORT_slasd4(aocl_int_t* n, aocl_int_t* i, float* d, float* z, float* delta, float* rho, float* sigma, float* work, aocl_int_t* info);
void LAPACK_EXPORT_slasd5(aocl_int_t* i, float* d, float* z, float* delta, float* rho, float* dsigma, float* work);
void LAPACK_EXPORT_slasd6(aocl_int_t* icompq, aocl_int_t* nl, aocl_int_t* nr, aocl_int_t* sqre, float* d, float* vf, float* vl, float* alpha, float* beta, aocl_int_t* idxq, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, aocl_int_t* ldgcol, float* givnum, aocl_int_t* ldgnum, float* poles, float* difl, float* difr, float* z, aocl_int_t* k, float* c, float* s, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_slasd7(aocl_int_t* icompq, aocl_int_t* nl, aocl_int_t* nr, aocl_int_t* sqre, aocl_int_t* k, float* d, float* z, float* zw, float* vf, float* vfw, float* vl, float* vlw, float* alpha, float* beta, float* dsigma, aocl_int_t* idx, aocl_int_t* idxp, aocl_int_t* idxq, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, aocl_int_t* ldgcol, float* givnum, aocl_int_t* ldgnum, float* c, float* s, aocl_int_t* info);
void LAPACK_EXPORT_slasd8(aocl_int_t* icompq, aocl_int_t* k, float* d, float* z, float* vf, float* vl, float* difl, float* difr, aocl_int_t* lddifr, float* dsigma, float* work, aocl_int_t* info);
void LAPACK_EXPORT_slasda(aocl_int_t* icompq, aocl_int_t* smlsiz, aocl_int_t* n, aocl_int_t* sqre, float* d, float* e, float* u, aocl_int_t* ldu, float* vt, aocl_int_t* k, float* difl, float* difr, float* z, float* poles, aocl_int_t* givptr, aocl_int_t* givcol, aocl_int_t* ldgcol, aocl_int_t* perm, float* givnum, float* c, float* s, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_slasdq(char* uplo, aocl_int_t* sqre, aocl_int_t* n, aocl_int_t* ncvt, aocl_int_t* nru, aocl_int_t* ncc, float* d, float* e, float* vt, aocl_int_t* ldvt, float* u, aocl_int_t* ldu, float* c, aocl_int_t* ldc, float* work, aocl_int_t* info);
void LAPACK_EXPORT_slasdt(aocl_int_t* n, aocl_int_t* lvl, aocl_int_t* nd, aocl_int_t* inode, aocl_int_t* ndiml, aocl_int_t* ndimr, aocl_int_t* msub);
void LAPACK_EXPORT_slaset(char* uplo, aocl_int_t* m, aocl_int_t* n, float* alpha, float* beta, float* a, aocl_int_t* lda);
void LAPACK_EXPORT_slasq1(aocl_int_t* n, float* d, float* e, float* work, aocl_int_t* info);
void LAPACK_EXPORT_slasq2(aocl_int_t* n, float* z, aocl_int_t* info);
void LAPACK_EXPORT_slasq3(aocl_int_t* i0, aocl_int_t* n0, float* z, aocl_int_t* pp, float* dmin, float* sigma, float* desig, float* qmax, aocl_int_t* nfail, aocl_int_t* iter, aocl_int_t* ndiv, logical* ieee, aocl_int_t* ttype, float* dmin1, float* dmin2, float* dn, float* dn1, float* dn2, float* g, float* tau);
void LAPACK_EXPORT_slasq4(aocl_int_t* i0, aocl_int_t* n0, float* z, aocl_int_t* pp, aocl_int_t* n0in, float* dmin, float* dmin1, float* dmin2, float* dn, float* dn1, float* dn2, float* tau, aocl_int_t* ttype, float* g);
void LAPACK_EXPORT_slasq5(aocl_int_t* i0, aocl_int_t* n0, float* z, aocl_int_t* pp, float* tau, float* sigma, float* dmin, float* dmin1, float* dmin2, float* dn, float* dnm1, float* dnm2, logical* ieee, float* eps);
void LAPACK_EXPORT_slasq6(aocl_int_t* i0, aocl_int_t* n0, float* z, aocl_int_t* pp, float* dmin, float* dmin1, float* dmin2, float* dn, float* dnm1, float* dnm2);
void LAPACK_EXPORT_slasr(char* side, char* pivot, char* direct, aocl_int_t* m, aocl_int_t* n, float* c, float* s, float* a, aocl_int_t* lda);
void LAPACK_EXPORT_slasrt(char* id, aocl_int_t* n, float* d, aocl_int_t* info);
void LAPACK_EXPORT_slassq(aocl_int_t* n, float* x, aocl_int_t* incx, float* scale, float* sumsq);
void LAPACK_EXPORT_slasv2(float* f, float* g, float* h, float* ssmin, float* ssmax, float* snr, float* csr, float* snl, float* csl);
void LAPACK_EXPORT_slaswlq(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, aocl_int_t* nb, float* a, aocl_int_t* lda, float* t, aocl_int_t* ldt, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_slaswp(aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* k1, aocl_int_t* k2, aocl_int_t* ipiv, aocl_int_t* incx);
void LAPACK_EXPORT_slasy2(logical* ltranl, logical* ltranr, aocl_int_t* isgn, aocl_int_t* n1, aocl_int_t* n2, float* tl, aocl_int_t* ldtl, float* tr, aocl_int_t* ldtr, float* b, aocl_int_t* ldb, float* scale, float* x, aocl_int_t* ldx, float* xnorm, aocl_int_t* info);
void LAPACK_EXPORT_sla_syamv(aocl_int_t* uplo, aocl_int_t* n, float* alpha, float* a, aocl_int_t* lda, float* x, aocl_int_t* incx, float* beta, float* y, aocl_int_t* incy);
void LAPACK_EXPORT_slasyf_aa(char* uplo, aocl_int_t* j1, aocl_int_t* m, aocl_int_t* nb, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* h, aocl_int_t* ldh, float* work);
void LAPACK_EXPORT_slasyf(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_slasyf_rk(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, float* a, aocl_int_t* lda, float* e, aocl_int_t* ipiv, float* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_slasyf_rook(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_sla_syrfsx_extended(aocl_int_t* prec_type, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, aocl_int_t* ipiv, logical* colequ, float* c, float* b, aocl_int_t* ldb, float* y, aocl_int_t* ldy, float* berr_out, aocl_int_t* n_norms, float* err_bnds_norm, float* err_bnds_comp, float* res, float* ayb, float* dy, float* y_tail, float* rcond, aocl_int_t* ithresh, float* rthresh, float* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_slatbs(char* uplo, char* trans, char* diag, char* normin, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, float* x, float* scale, float* cnorm, aocl_int_t* info);
void LAPACK_EXPORT_slatdf(aocl_int_t* ijob, aocl_int_t* n, float* z, aocl_int_t* ldz, float* rhs, float* rdsum, float* rdscal, aocl_int_t* ipiv, aocl_int_t* jpiv);
void LAPACK_EXPORT_slatps(char* uplo, char* trans, char* diag, char* normin, aocl_int_t* n, float* ap, float* x, float* scale, float* cnorm, aocl_int_t* info);
void LAPACK_EXPORT_slatrd(char* uplo, aocl_int_t* n, aocl_int_t* nb, float* a, aocl_int_t* lda, float* e, float* tau, float* w, aocl_int_t* ldw);
void LAPACK_EXPORT_slatrs(char* uplo, char* trans, char* diag, char* normin, aocl_int_t* n, float* a, aocl_int_t* lda, float* x, float* scale, float* cnorm, aocl_int_t* info);
void LAPACK_EXPORT_slatrz(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, float* a, aocl_int_t* lda, float* tau, float* work);
void LAPACK_EXPORT_slatsqr(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, aocl_int_t* nb, float* a, aocl_int_t* lda, float* t, aocl_int_t* ldt, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_slatzm(char* side, aocl_int_t* m, aocl_int_t* n, float* v, aocl_int_t* incv, float* tau, float* c1, float* c2, aocl_int_t* ldc, float* work);
void LAPACK_EXPORT_slauu2(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_slauum(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_sla_wwaddw(aocl_int_t* n, float* x, float* y, float* w);
void LAPACK_EXPORT_sopmtr(char* side, char* uplo, char* trans, aocl_int_t* m, aocl_int_t* n, float* ap, float* tau, float* c, aocl_int_t* ldc, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sorbdb1(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, float* x11, aocl_int_t* ldx11, float* x21, aocl_int_t* ldx21, float* theta, float* phi, float* taup1, float* taup2, float* tauq1, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorbdb2(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, float* x11, aocl_int_t* ldx11, float* x21, aocl_int_t* ldx21, float* theta, float* phi, float* taup1, float* taup2, float* tauq1, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorbdb3(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, float* x11, aocl_int_t* ldx11, float* x21, aocl_int_t* ldx21, float* theta, float* phi, float* taup1, float* taup2, float* tauq1, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorbdb4(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, float* x11, aocl_int_t* ldx11, float* x21, aocl_int_t* ldx21, float* theta, float* phi, float* taup1, float* taup2, float* tauq1, float* phantom, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorbdb5(aocl_int_t* m1, aocl_int_t* m2, aocl_int_t* n, float* x1, aocl_int_t* incx1, float* x2, aocl_int_t* incx2, float* q1, aocl_int_t* ldq1, float* q2, aocl_int_t* ldq2, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorbdb6(aocl_int_t* m1, aocl_int_t* m2, aocl_int_t* n, float* x1, aocl_int_t* incx1, float* x2, aocl_int_t* incx2, float* q1, aocl_int_t* ldq1, float* q2, aocl_int_t* ldq2, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorbdb(char* trans, char* signs, aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, float* x11, aocl_int_t* ldx11, float* x12, aocl_int_t* ldx12, float* x21, aocl_int_t* ldx21, float* x22, aocl_int_t* ldx22, float* theta, float* phi, float* taup1, float* taup2, float* tauq1, float* tauq2, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorg2l(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sorgl2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sorg2r(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sorgbr(char* vect, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorglq(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorgql(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorgqr(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorgr2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sorgrq(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorgtr(char* uplo, aocl_int_t* m, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorgtsqr(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, aocl_int_t* nb, float* a, aocl_int_t* lda, float* t, aocl_int_t* ldt, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorhr_col(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nb, float* a, aocl_int_t* lda, float* t, aocl_int_t* ldt, float* d, aocl_int_t* info);
void LAPACK_EXPORT_sorm22(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* n1, aocl_int_t* n2, float* q, aocl_int_t* ldq, float* c, aocl_int_t* ldc, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorm2l(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* c, aocl_int_t* ldc, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sorm2r(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* c, aocl_int_t* ldc, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sormbr(char* vect, char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* c, aocl_int_t* ldc, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sorml2(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* c, aocl_int_t* ldc, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sormlq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* c, aocl_int_t* ldc, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sormql(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* c, aocl_int_t* ldc, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sormqr(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* c, aocl_int_t* ldc, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sormr2(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* c, aocl_int_t* ldc, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sormr3(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, float* a, aocl_int_t* lda, float* tau, float* c, aocl_int_t* ldc, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sormrq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, float* a, aocl_int_t* lda, float* tau, float* c, aocl_int_t* ldc, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sormrz(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, float* a, aocl_int_t* lda, float* tau, float* c, aocl_int_t* ldc, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_sormtr(char* side, char* uplo, char* trans, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* tau, float* c, aocl_int_t* ldc, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_spbcon(char* uplo, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, float* anorm, float* rcond, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_spbequ(char* uplo, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, float* s, float* scond, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_spbrfs(char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, float* ab, aocl_int_t* ldab, float* afb, aocl_int_t* ldafb, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_spbstf(char* uplo, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, aocl_int_t* info);
void LAPACK_EXPORT_spbsv(char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, float* ab, aocl_int_t* ldab, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_spbsvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, float* ab, aocl_int_t* ldab, float* afb, aocl_int_t* ldafb, char* equed, float* s, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_spbtf2(char* uplo, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, aocl_int_t* info);
void LAPACK_EXPORT_spbtrf(char* uplo, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, aocl_int_t* info);
void LAPACK_EXPORT_spbtrs(char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, float* ab, aocl_int_t* ldab, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_spftrf(char* transr, char* uplo, aocl_int_t* n, float* a, aocl_int_t* info);
void LAPACK_EXPORT_spftri(char* transr, char* uplo, aocl_int_t* n, float* a, aocl_int_t* info);
void LAPACK_EXPORT_spftrs(char* transr, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_spocon(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* anorm, float* rcond, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_spoequb(aocl_int_t* n, float* a, aocl_int_t* lda, float* s, float* scond, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_spoequ(aocl_int_t* n, float* a, aocl_int_t* lda, float* s, float* scond, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_sporfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sporfsx(char* uplo, char* equed, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, float* s, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sposv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_sposvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, char* equed, float* s, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sposvxx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, char* equed, float* s, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* rpvgrw, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_spotf2(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_spotrf2(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_spotrf(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_spotri(char* uplo, aocl_int_t* n, float* buff_A, aocl_int_t* ldim_A, aocl_int_t* info);
void LAPACK_EXPORT_spotrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_sppcon(char* uplo, aocl_int_t* n, float* ap, float* anorm, float* rcond, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sppequ(char* uplo, aocl_int_t* n, float* ap, float* s, float* scond, float* amax, aocl_int_t* info);
void LAPACK_EXPORT_spprfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* ap, float* afp, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sppsv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* ap, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_sppsvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* ap, float* afp, char* equed, float* s, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_spptrf(char* uplo, aocl_int_t* n, float* ap, aocl_int_t* info);
void LAPACK_EXPORT_spptri(char* uplo, aocl_int_t* n, float* ap, aocl_int_t* info);
void LAPACK_EXPORT_spptrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* ap, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_spstf2(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* piv, aocl_int_t* rank, float* tol, float* work, aocl_int_t* info);
void LAPACK_EXPORT_spstrf(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* piv, aocl_int_t* rank, float* tol, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sptcon(aocl_int_t* n, float* d, float* e, float* anorm, float* rcond, float* work, aocl_int_t* info);
void LAPACK_EXPORT_spteqr(char* compz, aocl_int_t* n, float* d, float* e, float* z, aocl_int_t* ldz, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sptrfs(aocl_int_t* n, aocl_int_t* nrhs, float* d, float* e, float* df, float* ef, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* ferr, float* berr, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sptsv(aocl_int_t* n, aocl_int_t* nrhs, float* d, float* e, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_sptsvx(char* fact, aocl_int_t* n, aocl_int_t* nrhs, float* d, float* e, float* df, float* ef, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, float* work, aocl_int_t* info);
void LAPACK_EXPORT_spttrf(aocl_int_t* n, float* d, float* e, aocl_int_t* info);
void LAPACK_EXPORT_spttrs(aocl_int_t* n, aocl_int_t* nrhs, float* d, float* e, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_sptts2(aocl_int_t* n, aocl_int_t* nrhs, float* d, float* e, float* b, aocl_int_t* ldb);
void LAPACK_EXPORT_srscl(aocl_int_t* n, float* sa, float* sx, aocl_int_t* incx);
void LAPACK_EXPORT_ssb2st_kernels(char* uplo, logical* wantz, aocl_int_t* ttype, aocl_int_t* st, aocl_int_t* ed, aocl_int_t* sweep, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* ib, float* a, aocl_int_t* lda, float* v, float* tau, aocl_int_t* ldvt, float* work);
void LAPACK_EXPORT_ssbev_2stage(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssbev(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* info);
void LAPACK_EXPORT_ssbevd_2stage(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_ssbevd(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_ssbevx_2stage(char* jobz, char* range, char* uplo, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, float* q, aocl_int_t* ldq, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_ssbevx(char* jobz, char* range, char* uplo, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, float* q, aocl_int_t* ldq, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_ssbgst(char* vect, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, float* ab, aocl_int_t* ldab, float* bb, aocl_int_t* ldbb, float* x, aocl_int_t* ldx, float* work, aocl_int_t* info);
void LAPACK_EXPORT_ssbgv(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, float* ab, aocl_int_t* ldab, float* bb, aocl_int_t* ldbb, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* info);
void LAPACK_EXPORT_ssbgvd(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, float* ab, aocl_int_t* ldab, float* bb, aocl_int_t* ldbb, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_ssbgvx(char* jobz, char* range, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, float* ab, aocl_int_t* ldab, float* bb, aocl_int_t* ldbb, float* q, aocl_int_t* ldq, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_ssbtrd(char* vect, char* uplo, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, float* d, float* e, float* q, aocl_int_t* ldq, float* work, aocl_int_t* info);
void LAPACK_EXPORT_ssfrk(char* transr, char* uplo, char* trans, aocl_int_t* n, aocl_int_t* k, float* alpha, float* a, aocl_int_t* lda, float* beta, float* c__);
void LAPACK_EXPORT_sspcon(char* uplo, aocl_int_t* n, float* ap, aocl_int_t* ipiv, float* anorm, float* rcond, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sspev(char* jobz, char* uplo, aocl_int_t* n, float* ap, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sspevd(char* jobz, char* uplo, aocl_int_t* n, float* ap, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_sspevx(char* jobz, char* range, char* uplo, aocl_int_t* n, float* ap, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_sspgst(aocl_int_t* itype, char* uplo, aocl_int_t* n, float* ap, float* bp, aocl_int_t* info);
void LAPACK_EXPORT_sspgvd(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, float* ap, float* bp, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_sspgv(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, float* ap, float* bp, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sspgvx(aocl_int_t* itype, char* jobz, char* range, char* uplo, aocl_int_t* n, float* ap, float* bp, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_ssprfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* ap, float* afp, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sspsv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* ap, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_sspsvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* ap, float* afp, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_ssptrd(char* uplo, aocl_int_t* n, float* ap, float* d, float* e, float* tau, aocl_int_t* info);
void LAPACK_EXPORT_ssptrf(char* uplo, aocl_int_t* n, float* ap, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_ssptri(char* uplo, aocl_int_t* n, float* ap, aocl_int_t* ipiv, float* work, aocl_int_t* info);
void LAPACK_EXPORT_ssptrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* ap, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_sstebz(char* range, char* order, aocl_int_t* n, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, float* d, float* e, aocl_int_t* m, aocl_int_t* nsplit, float* w, aocl_int_t* iblock, aocl_int_t* isplit, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_sstedc(char* compz, aocl_int_t* n, float* d, float* e, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_sstegr(char* jobz, char* range, aocl_int_t* n, float* d, float* e, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, float* z, aocl_int_t* ldz, aocl_int_t* isuppz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_sstein(aocl_int_t* n, float* d, float* e, aocl_int_t* m, float* w, aocl_int_t* iblock, aocl_int_t* isplit, float* z, aocl_int_t* ldz, float* work, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_sstemr(char* jobz, char* range, aocl_int_t* n, float* d, float* e, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, aocl_int_t* m, float* w, float* z, aocl_int_t* ldz, aocl_int_t* nzc, aocl_int_t* isuppz, logical* tryrac, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_ssteqr(char* compz, aocl_int_t* n, float* d, float* e, float* z, aocl_int_t* ldz, float* work, aocl_int_t* info);
void LAPACK_EXPORT_ssteqr(char* jobz, aocl_int_t* n, float* d, float* e, float* z, aocl_int_t* ldz, float* work, aocl_int_t* info);
void LAPACK_EXPORT_ssterf(aocl_int_t* n, float* d, float* e, aocl_int_t* info);
void LAPACK_EXPORT_sstev(char* jobz, aocl_int_t* n, float* d, float* e, float* z, aocl_int_t* ldz, float* work, aocl_int_t* info);
void LAPACK_EXPORT_sstevd(char* jobz, aocl_int_t* n, float* d, float* e, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_sstevr(char* jobz, char* range, aocl_int_t* n, float* d, float* e, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, float* z, aocl_int_t* ldz, aocl_int_t* isuppz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_sstevx(char* jobz, char* range, aocl_int_t* n, float* d, float* e, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_ssycon_3(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* e, aocl_int_t* ipiv, float* anorm, float* rcond, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_ssycon(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* anorm, float* rcond, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_ssycon_rook(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* anorm, float* rcond, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_ssyconv(char* uplo, char* way, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* work, aocl_int_t* info);
void LAPACK_EXPORT_ssyconvf(char* uplo, char* way, aocl_int_t* n, float* a, aocl_int_t* lda, float* e, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_ssyconvf_rook(char* uplo, char* way, aocl_int_t* n, float* a, aocl_int_t* lda, float* e, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_ssyequb(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* s, float* scond, float* amax, float* work, aocl_int_t* info);
void LAPACK_EXPORT_ssyev_2stage(char* jobz, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* w, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssyev(char* jobz, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* w, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssyevd_2stage(char* jobz, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* w, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_ssyevd(char* jobz, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* w, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_ssyevr_2stage(char* jobz, char* range, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, float* z, aocl_int_t* ldz, aocl_int_t* isuppz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_ssyevr(char* jobz, char* range, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, float* z, aocl_int_t* ldz, aocl_int_t* isuppz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_ssyevx_2stage(char* jobz, char* range, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_ssyevx(char* jobz, char* range, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_ssygs2(aocl_int_t* itype, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_ssygst(aocl_int_t* itype, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_ssygv_2stage(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* w, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssygvd(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* w, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_ssygv(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* w, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssygvx(aocl_int_t* itype, char* jobz, char* range, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* vl, float* vu, aocl_int_t* il, aocl_int_t* iu, float* abstol, aocl_int_t* m, float* w, float* z, aocl_int_t* ldz, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_ssyrfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_ssyrfsx(char* uplo, char* equed, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, aocl_int_t* ipiv, float* s, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_ssysv_aa_2stage(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, float* b, aocl_int_t* ldb, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssysv_aa(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssysv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssysv_rk(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* e, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssysv_rook(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssysvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* ferr, float* berr, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_ssysvxx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* af, aocl_int_t* ldaf, aocl_int_t* ipiv, char* equed, float* s, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* rcond, float* rpvgrw, float* berr, aocl_int_t* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, aocl_int_t* nparams, float* params, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_ssyswapr(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* i1, aocl_int_t* i2);
void LAPACK_EXPORT_ssytd2(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* d, float* e, float* tau, aocl_int_t* info);
void LAPACK_EXPORT_ssytf2(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_ssytf2_rk(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* e, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_ssytf2_rook(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_ssytrd_2stage(char* vect, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* d, float* e, float* tau, float* hous2, aocl_int_t* lhous2, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssytrd(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* d, float* e, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssytrd_sb2st(char* stage1, char* vect, char* uplo, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, float* d, float* e, float* hous, aocl_int_t* lhous, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssytrd_sy2sb(char* uplo, aocl_int_t* n, aocl_int_t* kd, float* a, aocl_int_t* lda, float* ab, aocl_int_t* ldab, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssytrf_aa_2stage(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssytrf_aa(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssytrf(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssytrf_rk(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* e, aocl_int_t* ipiv, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssytrf_rook(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssytri2(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssytri2x(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* work, aocl_int_t* nb, aocl_int_t* info);
void LAPACK_EXPORT_ssytri_3(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* e, aocl_int_t* ipiv, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssytri_3x(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* e, aocl_int_t* ipiv, float* work, aocl_int_t* nb, aocl_int_t* info);
void LAPACK_EXPORT_ssytri(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* work, aocl_int_t* info);
void LAPACK_EXPORT_ssytri_rook(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* work, aocl_int_t* info);
void LAPACK_EXPORT_ssytrs2(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, float* work, aocl_int_t* info);
void LAPACK_EXPORT_ssytrs_3(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* e, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_ssytrs_aa_2stage(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_ssytrs_aa(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ssytrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_ssytrs_rook(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, aocl_int_t* ipiv, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_stbcon(char* norm, char* uplo, char* diag, aocl_int_t* n, aocl_int_t* kd, float* ab, aocl_int_t* ldab, float* rcond, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_stbrfs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, float* ab, aocl_int_t* ldab, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_stbtrs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, float* ab, aocl_int_t* ldab, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_stfsm(char* transr, char* side, char* uplo, char* trans, char* diag, aocl_int_t* m, aocl_int_t* n, float* alpha, float* a, float* b, aocl_int_t* ldb);
void LAPACK_EXPORT_stftri(char* transr, char* uplo, char* diag, aocl_int_t* n, float* a, aocl_int_t* info);
void LAPACK_EXPORT_stfttp(char* transr, char* uplo, aocl_int_t* n, float* arf, float* ap, aocl_int_t* info);
void LAPACK_EXPORT_stfttr(char* transr, char* uplo, aocl_int_t* n, float* arf, float* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_stgevc(char* side, char* howmny, logical* select, aocl_int_t* n, float* s, aocl_int_t* lds, float* p, aocl_int_t* ldp, float* vl, aocl_int_t* ldvl, float* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, float* work, aocl_int_t* info);
void LAPACK_EXPORT_stgex2(logical* wantq, logical* wantz, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* q, aocl_int_t* ldq, float* z, aocl_int_t* ldz, aocl_int_t* j1, aocl_int_t* n1, aocl_int_t* n2, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_stgexc(logical* wantq, logical* wantz, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* q, aocl_int_t* ldq, float* z, aocl_int_t* ldz, aocl_int_t* ifst, aocl_int_t* ilst, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_stgsen(aocl_int_t* ijob, logical* wantq, logical* wantz, logical* select, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* alphar, float* alphai, float* beta, float* q, aocl_int_t* ldq, float* z, aocl_int_t* ldz, aocl_int_t* m, float* pl, float* pr, float* dif, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_stgsja(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* tola, float* tolb, float* alpha, float* beta, float* u, aocl_int_t* ldu, float* v, aocl_int_t* ldv, float* q, aocl_int_t* ldq, float* work, aocl_int_t* ncycle, aocl_int_t* info);
void LAPACK_EXPORT_stgsna(char* job, char* howmny, logical* select, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* vl, aocl_int_t* ldvl, float* vr, aocl_int_t* ldvr, float* s, float* dif, aocl_int_t* mm, aocl_int_t* m, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_stgsy2(char* trans, aocl_int_t* ijob, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* c, aocl_int_t* ldc, float* d, aocl_int_t* ldd, float* e, aocl_int_t* lde, float* f, aocl_int_t* ldf, float* scale, float* rdsum, float* rdscal, aocl_int_t* iwork, aocl_int_t* pq, aocl_int_t* info);
void LAPACK_EXPORT_stgsyl(char* trans, aocl_int_t* ijob, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* c, aocl_int_t* ldc, float* d, aocl_int_t* ldd, float* e, aocl_int_t* lde, float* f, aocl_int_t* ldf, float* scale, float* dif, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_stpcon(char* norm, char* uplo, char* diag, aocl_int_t* n, float* ap, float* rcond, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_stplqt2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_stplqt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, aocl_int_t* mb, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* t, aocl_int_t* ldt, float* work, aocl_int_t* info);
void LAPACK_EXPORT_stpmlqt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, aocl_int_t* mb, float* v, aocl_int_t* ldv, float* t, aocl_int_t* ldt, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* work, aocl_int_t* info);
void LAPACK_EXPORT_stpmqrt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, aocl_int_t* nb, float* v, aocl_int_t* ldv, float* t, aocl_int_t* ldt, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* work, aocl_int_t* info);
void LAPACK_EXPORT_stpqrt2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_stpqrt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, aocl_int_t* nb, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* t, aocl_int_t* ldt, float* work, aocl_int_t* info);
void LAPACK_EXPORT_stprfb(char* side, char* trans, char* direct, char* storev, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, float* v, aocl_int_t* ldv, float* t, aocl_int_t* ldt, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* work, aocl_int_t* ldwork);
void LAPACK_EXPORT_stprfs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, float* ap, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_stptri(char* uplo, char* diag, aocl_int_t* n, float* ap, aocl_int_t* info);
void LAPACK_EXPORT_stptrs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, float* ap, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_stpttf(char* transr, char* uplo, aocl_int_t* n, float* ap, float* arf, aocl_int_t* info);
void LAPACK_EXPORT_stpttr(char* uplo, aocl_int_t* n, float* ap, float* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_strcon(char* norm, char* uplo, char* diag, aocl_int_t* n, float* a, aocl_int_t* lda, float* rcond, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_strevc3(char* side, char* howmny, logical* select, aocl_int_t* n, float* t, aocl_int_t* ldt, float* vl, aocl_int_t* ldvl, float* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_strevc(char* side, char* howmny, logical* select, aocl_int_t* n, float* t, aocl_int_t* ldt, float* vl, aocl_int_t* ldvl, float* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, float* work, aocl_int_t* info);
void LAPACK_EXPORT_strexc(char* compq, aocl_int_t* n, float* t, aocl_int_t* ldt, float* q, aocl_int_t* ldq, aocl_int_t* ifst, aocl_int_t* ilst, float* work, aocl_int_t* info);
void LAPACK_EXPORT_strrfs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* x, aocl_int_t* ldx, float* ferr, float* berr, float* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_strsen(char* job, char* compq, logical* select, aocl_int_t* n, float* t, aocl_int_t* ldt, float* q, aocl_int_t* ldq, float* wr, float* wi, aocl_int_t* m, float* s, float* sep, float* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_strsna(char* job, char* howmny, logical* select, aocl_int_t* n, float* t, aocl_int_t* ldt, float* vl, aocl_int_t* ldvl, float* vr, aocl_int_t* ldvr, float* s, float* sep, aocl_int_t* mm, aocl_int_t* m, float* work, aocl_int_t* ldwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_strsyl(char* trana, char* tranb, aocl_int_t* isgn, aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, float* c, aocl_int_t* ldc, float* scale, aocl_int_t* info);
void LAPACK_EXPORT_strti2(char* uplo, char* diag, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_strtri(char* uplo, char* diag, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_strtrs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, float* a, aocl_int_t* lda, float* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_strttf(char* transr, char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* arf, aocl_int_t* info);
void LAPACK_EXPORT_strttp(char* uplo, aocl_int_t* n, float* a, aocl_int_t* lda, float* ap, aocl_int_t* info);
void LAPACK_EXPORT_stzrqf(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* tau, aocl_int_t* info);
void LAPACK_EXPORT_stzrzf(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, float* tau, float* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zbbcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, double* theta, double* phi, dcomplex* u1, aocl_int_t* ldu1, dcomplex* u2, aocl_int_t* ldu2, dcomplex* v1t, aocl_int_t* ldv1t, dcomplex* v2t, aocl_int_t* ldv2t, double* b11d, double* b11e, double* b12d, double* b12e, double* b21d, double* b21e, double* b22d, double* b22e, double* rwork, aocl_int_t* lrwork, aocl_int_t* info);
void LAPACK_EXPORT_zbdsqr(char* uplo, aocl_int_t* n, aocl_int_t* ncvt, aocl_int_t* nru, aocl_int_t* ncc, double* d, double* e, dcomplex* vt, aocl_int_t* ldvt, dcomplex* u, aocl_int_t* ldu, dcomplex* c, aocl_int_t* ldc, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zcgesv(aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, dcomplex* work, scomplex* swork, double* rwork, aocl_int_t* iter, aocl_int_t* info);
void LAPACK_EXPORT_zcposv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, dcomplex* work, scomplex* swork, double* rwork, aocl_int_t* iter, aocl_int_t* info);
void LAPACK_EXPORT_zdrscl(aocl_int_t* n, double* sa, dcomplex* sx, aocl_int_t* incx);
void LAPACK_EXPORT_zgbbrd(char* vect, aocl_int_t* m, aocl_int_t* n, aocl_int_t* ncc, aocl_int_t* kl, aocl_int_t* ku, dcomplex* ab, aocl_int_t* ldab, double* d, double* e, dcomplex* q, aocl_int_t* ldq, dcomplex* pt, aocl_int_t* ldpt, dcomplex* c, aocl_int_t* ldc, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgbcon(char* norm, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, dcomplex* ab, aocl_int_t* ldab, aocl_int_t* ipiv, double* anorm, double* rcond, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgbequb(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, dcomplex* ab, aocl_int_t* ldab, double* r, double* c, double* rowcnd, double* colcnd, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_zgbequ(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, dcomplex* ab, aocl_int_t* ldab, double* r, double* c, double* rowcnd, double* colcnd, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_zgbrfs(char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, dcomplex* ab, aocl_int_t* ldab, dcomplex* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgbrfsx(char* trans, char* equed, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, dcomplex* ab, aocl_int_t* ldab, dcomplex* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, double* r, double* c, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgbsv(aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, dcomplex* ab, aocl_int_t* ldab, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zgbsvx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, dcomplex* ab, aocl_int_t* ldab, dcomplex* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, char* equed, double* r, double* c, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgbsvxx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, dcomplex* ab, aocl_int_t* ldab, dcomplex* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, char* equed, double* r, double* c, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* rpvgrw, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgbtf2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, dcomplex* ab, aocl_int_t* ldab, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zgbtrf(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, dcomplex* ab, aocl_int_t* ldab, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zgbtrs(char* trans, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, dcomplex* ab, aocl_int_t* ldab, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zgebak(char* job, char* side, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, double* scale, aocl_int_t* m, dcomplex* v, aocl_int_t* ldv, aocl_int_t* info);
void LAPACK_EXPORT_zgebal(char* job, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ilo, aocl_int_t* ihi, double* scale, aocl_int_t* info);
void LAPACK_EXPORT_zgebd2(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* d, double* e, dcomplex* tauq, dcomplex* taup, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zgebrd(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* d, double* e, dcomplex* tauq, dcomplex* taup, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgecon(char* norm, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* anorm, double* rcond, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgeequb(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_zgeequ(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_zgees(char* jobvs, char* sort, L_fp select, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* sdim, dcomplex* w, dcomplex* vs, aocl_int_t* ldvs, dcomplex* work, aocl_int_t* lwork, double* rwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_zgeesx(char* jobvs, char* sort, L_fp select, char* sense, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* sdim, dcomplex* w, dcomplex* vs, aocl_int_t* ldvs, double* rconde, double* rcondv, dcomplex* work, aocl_int_t* lwork, double* rwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_zgeev(char* jobvl, char* jobvr, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* w, dcomplex* vl, aocl_int_t* ldvl, dcomplex* vr, aocl_int_t* ldvr, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgeevx(char* balanc, char* jobvl, char* jobvr, char* sense, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* w, dcomplex* vl, aocl_int_t* ldvl, dcomplex* vr, aocl_int_t* ldvr, aocl_int_t* ilo, aocl_int_t* ihi, double* scale, double* abnrm, double* rconde, double* rcondv, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgegs(char* jobvsl, char* jobvsr, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vsl, aocl_int_t* ldvsl, dcomplex* vsr, aocl_int_t* ldvsr, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgegv(char* jobvl, char* jobvr, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vl, aocl_int_t* ldvl, dcomplex* vr, aocl_int_t* ldvr, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgehd2(aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zgehrd(aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgejsv(char* joba, char* jobu, char* jobv, char* jobr, char* jobt, char* jobp, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* sva, dcomplex* u, aocl_int_t* ldu, dcomplex* v, aocl_int_t* ldv, dcomplex* cwork, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_zgelq2(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zgelqf(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgelq(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* tsize, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgelqt3(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_zgelqt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* ldt, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zgels(char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgelsd(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, double* s, double* rcond, aocl_int_t* rank, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_zgelss(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, double* s, double* rcond, aocl_int_t* rank, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgelsx(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, aocl_int_t* jpvt, double* rcond, aocl_int_t* rank, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgelsy(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, aocl_int_t* jpvt, double* rcond, aocl_int_t* rank, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgemlq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* tsize, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgemlqt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* mb, dcomplex* v, aocl_int_t* ldv, dcomplex* t, aocl_int_t* ldt, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zgemqr(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* tsize, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgemqrt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* nb, dcomplex* v, aocl_int_t* ldv, dcomplex* t, aocl_int_t* ldt, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zgeql2(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zgeqlf(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgeqp3(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* jpvt, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgeqpf(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* jpvt, dcomplex* tau, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgeqr2(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zgeqr2p(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zgeqrf(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgeqrfp(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgeqr(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* tsize, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgeqrt2(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_zgeqrt3(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_zgeqrt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nb, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* ldt, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zgerfs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgerfsx(char* trans, char* equed, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, double* r, double* c, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgerq2(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zgerqf(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgesc2(aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* rhs, aocl_int_t* ipiv, aocl_int_t* jpiv, double* scale);
void LAPACK_EXPORT_zgesdd(char* jobz, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* s, dcomplex* u, aocl_int_t* ldu, dcomplex* vt, aocl_int_t* ldvt, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_zgesvd(char* jobu, char* jobv, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* s, dcomplex* u, aocl_int_t* ldu, dcomplex* vt, aocl_int_t* ldvt, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgesvdq(char* joba, char* jobp, char* jobr, char* jobu, char* jobv, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* s, dcomplex* u, aocl_int_t* ldu, dcomplex* v, aocl_int_t* ldv, aocl_int_t* numrank, aocl_int_t* iwork, aocl_int_t* liwork, dcomplex* cwork, aocl_int_t* lcwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* info);
void LAPACK_EXPORT_zgesvdx(char* jobu, char* jobvt, char* range, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, aocl_int_t* ns, double* s, dcomplex* u, aocl_int_t* ldu, dcomplex* vt, aocl_int_t* ldvt, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_zgesv(aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zgesvj(char* joba, char* jobu, char* jobv, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* sva, aocl_int_t* mv, dcomplex* v, aocl_int_t* ldv, dcomplex* cwork, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* info);
void LAPACK_EXPORT_zgesvx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, char* equed, double* r, double* c, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgesvxx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, char* equed, double* r, double* c, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* rpvgrw, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgetc2(aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* jpiv, aocl_int_t* info);
void LAPACK_EXPORT_zgetf2(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zgetrf2(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zgetrf(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zgetri(aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgetrs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zgetsls(char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zggbak(char* job, char* side, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, double* lscale, double* rscale, aocl_int_t* m, dcomplex* v, aocl_int_t* ldv, aocl_int_t* info);
void LAPACK_EXPORT_zggbal(char* job, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, aocl_int_t* ilo, aocl_int_t* ihi, double* lscale, double* rscale, double* work, aocl_int_t* info);
void LAPACK_EXPORT_zgges3(char* jobvsl, char* jobvsr, char* sort, L_fp selctg, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, aocl_int_t* sdim, dcomplex* alpha, dcomplex* beta, dcomplex* vsl, aocl_int_t* ldvsl, dcomplex* vsr, aocl_int_t* ldvsr, dcomplex* work, aocl_int_t* lwork, double* rwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_zgges(char* jobvsl, char* jobvsr, char* sort, L_fp selctg, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, aocl_int_t* sdim, dcomplex* alpha, dcomplex* beta, dcomplex* vsl, aocl_int_t* ldvsl, dcomplex* vsr, aocl_int_t* ldvsr, dcomplex* work, aocl_int_t* lwork, double* rwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_zggesx(char* jobvsl, char* jobvsr, char* sort, L_fp selctg, char* sense, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, aocl_int_t* sdim, dcomplex* alpha, dcomplex* beta, dcomplex* vsl, aocl_int_t* ldvsl, dcomplex* vsr, aocl_int_t* ldvsr, double* rconde, double* rcondv, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* iwork, aocl_int_t* liwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_zggev3(char* jobvl, char* jobvr, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vl, aocl_int_t* ldvl, dcomplex* vr, aocl_int_t* ldvr, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zggev(char* jobvl, char* jobvr, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vl, aocl_int_t* ldvl, dcomplex* vr, aocl_int_t* ldvr, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zggevx(char* balanc, char* jobvl, char* jobvr, char* sense, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vl, aocl_int_t* ldvl, dcomplex* vr, aocl_int_t* ldvr, aocl_int_t* ilo, aocl_int_t* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, double* rconde, double* rcondv, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* iwork, logical* bwork, aocl_int_t* info);
void LAPACK_EXPORT_zggglm(aocl_int_t* n, aocl_int_t* m, aocl_int_t* p, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* d, dcomplex* x, dcomplex* y, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgghd3(char* compq, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* q, aocl_int_t* ldq, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgghrd(char* compq, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* q, aocl_int_t* ldq, dcomplex* z, aocl_int_t* ldz, aocl_int_t* info);
void LAPACK_EXPORT_zgglse(aocl_int_t* m, aocl_int_t* n, aocl_int_t* p, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* c, dcomplex* d, dcomplex* x, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zggqrf(aocl_int_t* n, aocl_int_t* m, aocl_int_t* p, dcomplex* a, aocl_int_t* lda, dcomplex* taua, dcomplex* b, aocl_int_t* ldb, dcomplex* taub, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zggrqf(aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* taua, dcomplex* b, aocl_int_t* ldb, dcomplex* taub, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zggsvd3(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* n, aocl_int_t* p, aocl_int_t* k, aocl_int_t* l, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, double* alpha, double* beta, dcomplex* u, aocl_int_t* ldu, dcomplex* v, aocl_int_t* ldv, dcomplex* q, aocl_int_t* ldq, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_zggsvd(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* n, aocl_int_t* p, aocl_int_t* k, aocl_int_t* l, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, double* alpha, double* beta, dcomplex* u, aocl_int_t* ldu, dcomplex* v, aocl_int_t* ldv, dcomplex* q, aocl_int_t* ldq, dcomplex* work, double* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_zggsvp3(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, double* tola, double* tolb, aocl_int_t* k, aocl_int_t* l, dcomplex* u, aocl_int_t* ldu, dcomplex* v, aocl_int_t* ldv, dcomplex* q, aocl_int_t* ldq, aocl_int_t* iwork, double* rwork, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zggsvp(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, double* tola, double* tolb, aocl_int_t* k, aocl_int_t* l, dcomplex* u, aocl_int_t* ldu, dcomplex* v, aocl_int_t* ldv, dcomplex* q, aocl_int_t* ldq, aocl_int_t* iwork, double* rwork, dcomplex* tau, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zgsvj0(char* jobv, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* d, double* sva, aocl_int_t* mv, dcomplex* v, aocl_int_t* ldv, double* eps, double* sfmin, double* tol, aocl_int_t* nsweep, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgsvj1(char* jobv, aocl_int_t* m, aocl_int_t* n, aocl_int_t* n1, dcomplex* a, aocl_int_t* lda, dcomplex* d, double* sva, aocl_int_t* mv, dcomplex* v, aocl_int_t* ldv, double* eps, double* sfmin, double* tol, aocl_int_t* nsweep, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zgtcon(char* norm, aocl_int_t* n, dcomplex* dl, dcomplex* d, dcomplex* du, dcomplex* du2, aocl_int_t* ipiv, double* anorm, double* rcond, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zgtrfs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* dl, dcomplex* d, dcomplex* du, dcomplex* dlf, dcomplex* df, dcomplex* duf, dcomplex* du2, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgtsv(aocl_int_t* n, aocl_int_t* nrhs, dcomplex* dl, dcomplex* d, dcomplex* du, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zgtsvx(char* fact, char* trans, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* dl, dcomplex* d, dcomplex* du, dcomplex* dlf, dcomplex* df, dcomplex* duf, dcomplex* du2, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zgttrf(aocl_int_t* n, dcomplex* dl, dcomplex* d, dcomplex* du, dcomplex* du2, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zgttrs(char* trans, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* dl, dcomplex* d, dcomplex* du, dcomplex* du2, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zgtts2(aocl_int_t* itrans, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* dl, dcomplex* d, dcomplex* du, dcomplex* du2, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb);
void LAPACK_EXPORT_zhb2st_kernels(char* uplo, logical* wantz, aocl_int_t* ttype, aocl_int_t* st, aocl_int_t* ed, aocl_int_t* sweep, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* ib, dcomplex* a, aocl_int_t* lda, dcomplex* v, dcomplex* tau, aocl_int_t* ldvt, dcomplex* work);
void LAPACK_EXPORT_zhbev_2stage(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zhbev(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zhbevd_2stage(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_zhbevd(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_zhbevx_2stage(char* jobz, char* range, char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, dcomplex* q, aocl_int_t* ldq, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_zhbevx(char* jobz, char* range, char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, dcomplex* q, aocl_int_t* ldq, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, double* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_zhbgst(char* vect, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, dcomplex* ab, aocl_int_t* ldab, dcomplex* bb, aocl_int_t* ldbb, dcomplex* x, aocl_int_t* ldx, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zhbgv(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, dcomplex* ab, aocl_int_t* ldab, dcomplex* bb, aocl_int_t* ldbb, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zhbgvd(char* jobz, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, dcomplex* ab, aocl_int_t* ldab, dcomplex* bb, aocl_int_t* ldbb, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_zhbgvx(char* jobz, char* range, char* uplo, aocl_int_t* n, aocl_int_t* ka, aocl_int_t* kb, dcomplex* ab, aocl_int_t* ldab, dcomplex* bb, aocl_int_t* ldbb, dcomplex* q, aocl_int_t* ldq, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, double* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_zhbtrd(char* vect, char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, double* d, double* e, dcomplex* q, aocl_int_t* ldq, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zhecon_3(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, double* anorm, double* rcond, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zhecon(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, double* anorm, double* rcond, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zhecon_rook(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, double* anorm, double* rcond, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zheequb(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* s, double* scond, double* amax, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zheev_2stage(char* jobz, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* w, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zheev(char* jobz, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* w, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zheevd_2stage(char* jobz, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* w, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_zheevd(char* jobz, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* w, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_zheevr_2stage(char* jobz, char* range, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, dcomplex* z, aocl_int_t* ldz, aocl_int_t* isuppz, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_zheevr(char* jobz, char* range, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, dcomplex* z, aocl_int_t* ldz, aocl_int_t* isuppz, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_zheevx_2stage(char* jobz, char* range, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_zheevx(char* jobz, char* range, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_zhegs2(aocl_int_t* itype, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zhegst(aocl_int_t* itype, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zhegv_2stage(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, double* w, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zhegvd(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, double* w, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_zhegv(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, double* w, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zhegvx(aocl_int_t* itype, char* jobz, char* range, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_zherfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zherfsx(char* uplo, char* equed, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, double* s, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zhesv_aa_2stage(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhesv_aa(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhesv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhesv_rk(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhesv_rook(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhesvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zhesvxx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, char* equed, double* s, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* rpvgrw, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zheswapr(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* i1, aocl_int_t* i2);
void LAPACK_EXPORT_zhetd2(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* d, double* e, dcomplex* tau, aocl_int_t* info);
void LAPACK_EXPORT_zhetf2(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zhetf2_rk(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zhetf2_rook(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zhetrd_2stage(char* vect, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* d, double* e, dcomplex* tau, dcomplex* hous2, aocl_int_t* lhous2, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhetrd(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* d, double* e, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhetrd_hb2st(char* stage1, char* vect, char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, double* d, double* e, dcomplex* hous, aocl_int_t* lhous, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhetrd_he2hb(char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* a, aocl_int_t* lda, dcomplex* ab, aocl_int_t* ldab, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhetrf_aa_2stage(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhetrf_aa(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhetrf(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhetrf_rk(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhetrf_rook(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhetri2(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhetri2x(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* nb, aocl_int_t* info);
void LAPACK_EXPORT_zhetri_3(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhetri_3x(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* nb, aocl_int_t* info);
void LAPACK_EXPORT_zhetri(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zhetri_rook(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zhetrs2(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zhetrs_3(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zhetrs_aa_2stage(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zhetrs_aa(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zhetrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zhetrs_rook(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zhfrk(char* transr, char* uplo, char* trans, aocl_int_t* n, aocl_int_t* k, double* alpha, dcomplex* a, aocl_int_t* lda, double* beta, dcomplex* c__);
void LAPACK_EXPORT_zhgeqz(char* job, char* compq, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, dcomplex* h, aocl_int_t* ldh, dcomplex* t, aocl_int_t* ldt, dcomplex* alpha, dcomplex* beta, dcomplex* q, aocl_int_t* ldq, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zhpcon(char* uplo, aocl_int_t* n, dcomplex* ap, aocl_int_t* ipiv, double* anorm, double* rcond, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zhpev(char* jobz, char* uplo, aocl_int_t* n, dcomplex* ap, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zhpevd(char* jobz, char* uplo, aocl_int_t* n, dcomplex* ap, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_zhpevx(char* jobz, char* range, char* uplo, aocl_int_t* n, dcomplex* ap, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, double* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_zhpgst(aocl_int_t* itype, char* uplo, aocl_int_t* n, dcomplex* ap, dcomplex* bp, aocl_int_t* info);
void LAPACK_EXPORT_zhpgvd(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, dcomplex* ap, dcomplex* bp, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_zhpgv(aocl_int_t* itype, char* jobz, char* uplo, aocl_int_t* n, dcomplex* ap, dcomplex* bp, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zhpgvx(aocl_int_t* itype, char* jobz, char* range, char* uplo, aocl_int_t* n, dcomplex* ap, dcomplex* bp, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, double* rwork, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_zhprfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* ap, dcomplex* afp, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zhpsv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* ap, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zhpsvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* ap, dcomplex* afp, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zhptrd(char* uplo, aocl_int_t* n, dcomplex* ap, double* d, double* e, dcomplex* tau, aocl_int_t* info);
void LAPACK_EXPORT_zhptrf(char* uplo, aocl_int_t* n, dcomplex* ap, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zhptri(char* uplo, aocl_int_t* n, dcomplex* ap, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zhptrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* ap, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zhsein(char* side, char* eigsrc, char* initv, logical* select, aocl_int_t* n, dcomplex* h, aocl_int_t* ldh, dcomplex* w, dcomplex* vl, aocl_int_t* ldvl, dcomplex* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, dcomplex* work, double* rwork, aocl_int_t* ifaill, aocl_int_t* ifailr, aocl_int_t* info);
void LAPACK_EXPORT_zhseqr(char* job, char* compz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, dcomplex* h, aocl_int_t* ldh, dcomplex* w, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zlabrd(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nb, dcomplex* a, aocl_int_t* lda, double* d, double* e, dcomplex* tauq, dcomplex* taup, dcomplex* x, aocl_int_t* ldx, dcomplex* y, aocl_int_t* ldy);
void LAPACK_EXPORT_zlacgv(aocl_int_t* n, dcomplex* x, aocl_int_t* incx);
void LAPACK_EXPORT_zlacn2(aocl_int_t* n, dcomplex* v, dcomplex* x, double* est, aocl_int_t* kase, aocl_int_t* isave);
void LAPACK_EXPORT_zlacon(aocl_int_t* n, dcomplex* v, dcomplex* x, double* est, aocl_int_t* kase);
void LAPACK_EXPORT_zlacp2(char* uplo, aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb);
void LAPACK_EXPORT_zlacpy(char* uplo, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb);
void LAPACK_EXPORT_zlacrm(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* b, aocl_int_t* ldb, dcomplex* c, aocl_int_t* ldc, double* rwork);
void LAPACK_EXPORT_zlacrt(aocl_int_t* n, dcomplex* cx, aocl_int_t* incx, dcomplex* cy, aocl_int_t* incy, dcomplex* c, dcomplex* s);
void LAPACK_EXPORT_zlaed0(aocl_int_t* qsiz, aocl_int_t* n, double* d, double* e, dcomplex* q, aocl_int_t* ldq, dcomplex* qstore, aocl_int_t* ldqs, double* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_zlaed7(aocl_int_t* n, aocl_int_t* cutpnt, aocl_int_t* qsiz, aocl_int_t* tlvls, aocl_int_t* curlvl, aocl_int_t* curpbm, double* d, dcomplex* q, aocl_int_t* ldq, double* rho, aocl_int_t* indxq, double* qstore, aocl_int_t* qptr, aocl_int_t* prmptr, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, double* givnum, dcomplex* work, double* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_zlaed8(aocl_int_t* k, aocl_int_t* n, aocl_int_t* qsiz, dcomplex* q, aocl_int_t* ldq, double* d, double* rho, aocl_int_t* cutpnt, double* z, double* dlamda, dcomplex* q2, aocl_int_t* ldq2, double* w, aocl_int_t* indxp, aocl_int_t* indx, aocl_int_t* indxq, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, double* givnum, aocl_int_t* info);
void LAPACK_EXPORT_zlaein(logical* rightv, logical* noinit, aocl_int_t* n, dcomplex* h, aocl_int_t* ldh, dcomplex* w, dcomplex* v, dcomplex* b, aocl_int_t* ldb, double* rwork, double* eps3, double* smlnum, aocl_int_t* info);
void LAPACK_EXPORT_zlaesy(dcomplex* a, dcomplex* b, dcomplex* c, dcomplex* rt1, dcomplex* rt2, dcomplex* evscal, dcomplex* cs1, dcomplex* sn1);
void LAPACK_EXPORT_zlaev2(dcomplex* a, dcomplex* b, dcomplex* c, double* rt1, double* rt2, double* cs1, dcomplex* sn1);
void LAPACK_EXPORT_zlag2c(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, scomplex* sa, aocl_int_t* ldsa, aocl_int_t* info);
void LAPACK_EXPORT_zla_gbamv(aocl_int_t* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, double* alpha, dcomplex* ab, aocl_int_t* ldab, dcomplex* x, aocl_int_t* incx, double* beta, double* y, aocl_int_t* incy);
void LAPACK_EXPORT_zla_gbrfsx_extended(aocl_int_t* prec_type, aocl_int_t* trans_type, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, aocl_int_t* nrhs, dcomplex* ab, aocl_int_t* ldab, dcomplex* afb, aocl_int_t* ldafb, aocl_int_t* ipiv, logical* colequ, double* c, dcomplex* b, aocl_int_t* ldb, dcomplex* y, aocl_int_t* ldy, double* berr_out, aocl_int_t* n_norms, double* err_bnds_norm, double* err_bnds_comp, dcomplex* res, double* ayb, dcomplex* dy, dcomplex* y_tail, double* rcond, aocl_int_t* ithresh, double* rthresh, double* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_zla_geamv(aocl_int_t* trans, aocl_int_t* m, aocl_int_t* n, double* alpha, dcomplex* a, aocl_int_t* lda, dcomplex* x, aocl_int_t* incx, double* beta, double* y, aocl_int_t* incy);
void LAPACK_EXPORT_zla_gerfsx_extended(aocl_int_t* prec_type, aocl_int_t* trans_type, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, logical* colequ, double* c, dcomplex* b, aocl_int_t* ldb, dcomplex* y, aocl_int_t* ldy, double* berr_out, aocl_int_t* n_norms, double* errs_n, double* errs_c, dcomplex* res, double* ayb, dcomplex* dy, dcomplex* y_tail, double* rcond, aocl_int_t* ithresh, double* rthresh, double* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_zlags2(logical* upper, double* a1, dcomplex* a2, double* a3, double* b1, dcomplex* b2, double* b3, double* csu, dcomplex* snu, double* csv, dcomplex* snv, double* csq, dcomplex* snq);
void LAPACK_EXPORT_zlagtm(char* trans, aocl_int_t* n, aocl_int_t* nrhs, double* alpha, dcomplex* dl, dcomplex* d, dcomplex* du, dcomplex* x, aocl_int_t* ldx, double* beta, dcomplex* b, aocl_int_t* ldb);
void LAPACK_EXPORT_zla_heamv(aocl_int_t* uplo, aocl_int_t* n, double* alpha, dcomplex* a, aocl_int_t* lda, dcomplex* x, aocl_int_t* incx, double* beta, double* y, aocl_int_t* incy);
void LAPACK_EXPORT_zlahef_aa(char* uplo, aocl_int_t* j1, aocl_int_t* m, aocl_int_t* nb, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* h, aocl_int_t* ldh, dcomplex* work);
void LAPACK_EXPORT_zlahef(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_zlahef_rk(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, dcomplex* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_zlahef_rook(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_zla_herfsx_extended(aocl_int_t* prec_type, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, logical* colequ, double* c, dcomplex* b, aocl_int_t* ldb, dcomplex* y, aocl_int_t* ldy, double* berr_out, aocl_int_t* n_norms, double* err_bnds_norm, double* err_bnds_comp, dcomplex* res, double* ayb, dcomplex* dy, dcomplex* y_tail, double* rcond, aocl_int_t* ithresh, double* rthresh, double* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_zlahqr(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, dcomplex* h, aocl_int_t* ldh, dcomplex* w, aocl_int_t* iloz, aocl_int_t* ihiz, dcomplex* z, aocl_int_t* ldz, aocl_int_t* info);
void LAPACK_EXPORT_zlahr2(aocl_int_t* n, aocl_int_t* k, aocl_int_t* nb, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* t, aocl_int_t* ldt, dcomplex* y, aocl_int_t* ldy);
void LAPACK_EXPORT_zlahrd(aocl_int_t* n, aocl_int_t* k, aocl_int_t* nb, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* t, aocl_int_t* ldt, dcomplex* y, aocl_int_t* ldy);
void LAPACK_EXPORT_zlaic1(aocl_int_t* job, aocl_int_t* j, dcomplex* x, double* sest, dcomplex* w, dcomplex* gamma, double* sestpr, dcomplex* s, dcomplex* c__);
void LAPACK_EXPORT_zla_lin_berr(aocl_int_t* n, aocl_int_t* nz, aocl_int_t* nrhs, dcomplex* res, double* ayb, double* berr);
void LAPACK_EXPORT_zlals0(aocl_int_t* icompq, aocl_int_t* nl, aocl_int_t* nr, aocl_int_t* sqre, aocl_int_t* nrhs, dcomplex* b, aocl_int_t* ldb, dcomplex* bx, aocl_int_t* ldbx, aocl_int_t* perm, aocl_int_t* givptr, aocl_int_t* givcol, aocl_int_t* ldgcol, double* givnum, aocl_int_t* ldgnum, double* poles, double* difl, double* difr, double* z, aocl_int_t* k, double* c, double* s, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zlalsa(aocl_int_t* icompq, aocl_int_t* smlsiz, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* b, aocl_int_t* ldb, dcomplex* bx, aocl_int_t* ldbx, double* u, aocl_int_t* ldu, double* vt, aocl_int_t* k, double* difl, double* difr, double* z, double* poles, aocl_int_t* givptr, aocl_int_t* givcol, aocl_int_t* ldgcol, aocl_int_t* perm, double* givnum, double* c, double* s, double* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_zlalsd(char* uplo, aocl_int_t* smlsiz, aocl_int_t* n, aocl_int_t* nrhs, double* d, double* e, dcomplex* b, aocl_int_t* ldb, double* rcond, aocl_int_t* rank, dcomplex* work, double* rwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_zlamswlq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* mb, aocl_int_t* nb, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* ldt, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zlamtsqr(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* mb, aocl_int_t* nb, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* ldt, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zlapll(aocl_int_t* n, dcomplex* x, aocl_int_t* incx, dcomplex* y, aocl_int_t* incy, double* ssmin);
void LAPACK_EXPORT_zlapmr(logical* forwrd, aocl_int_t* m, aocl_int_t* n, dcomplex* x, aocl_int_t* ldx, aocl_int_t* k);
void LAPACK_EXPORT_zlapmt(logical* forwrd, aocl_int_t* m, aocl_int_t* n, dcomplex* x, aocl_int_t* ldx, aocl_int_t* k);
void LAPACK_EXPORT_zla_porfsx_extended(aocl_int_t* prec_type, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, logical* colequ, double* c, dcomplex* b, aocl_int_t* ldb, dcomplex* y, aocl_int_t* ldy, double* berr_out, aocl_int_t* n_norms, double* err_bnds_norm, double* err_bnds_comp, dcomplex* res, double* ayb, dcomplex* dy, dcomplex* y_tail, double* rcond, aocl_int_t* ithresh, double* rthresh, double* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_zlaqgb(aocl_int_t* m, aocl_int_t* n, aocl_int_t* kl, aocl_int_t* ku, dcomplex* ab, aocl_int_t* ldab, double* r, double* c, double* rowcnd, double* colcnd, double* amax, char* equed);
void LAPACK_EXPORT_zlaqge(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, char* equed);
void LAPACK_EXPORT_zlaqhb(char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, double* s, double* scond, double* amax, char* equed);
void LAPACK_EXPORT_zlaqhe(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* s, double* scond, double* amax, char* equed);
void LAPACK_EXPORT_zlaqhp(char* uplo, aocl_int_t* n, dcomplex* ap, double* s, double* scond, double* amax, char* equed);
void LAPACK_EXPORT_zlaqp2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* offset, dcomplex* a, aocl_int_t* lda, aocl_int_t* jpvt, dcomplex* tau, double* vn1, double* vn2, dcomplex* work);
void LAPACK_EXPORT_zlaqps(aocl_int_t* m, aocl_int_t* n, aocl_int_t* offset, aocl_int_t* nb, aocl_int_t* kb, dcomplex* a, aocl_int_t* lda, aocl_int_t* jpvt, dcomplex* tau, double* vn1, double* vn2, dcomplex* auxv, dcomplex* f, aocl_int_t* ldf);
void LAPACK_EXPORT_zlaqr0(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, dcomplex* h, aocl_int_t* ldh, dcomplex* w, aocl_int_t* iloz, aocl_int_t* ihiz, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zlaqr1(aocl_int_t* n, dcomplex* h, aocl_int_t* ldh, dcomplex* s1, dcomplex* s2, dcomplex* v);
void LAPACK_EXPORT_zlaqr2(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ktop, aocl_int_t* kbot, aocl_int_t* nw, dcomplex* h, aocl_int_t* ldh, aocl_int_t* iloz, aocl_int_t* ihiz, dcomplex* z, aocl_int_t* ldz, aocl_int_t* ns, aocl_int_t* nd, dcomplex* sh, dcomplex* v, aocl_int_t* ldv, aocl_int_t* nh, dcomplex* t, aocl_int_t* ldt, aocl_int_t* nv, dcomplex* wv, aocl_int_t* ldwv, dcomplex* work, aocl_int_t* lwork);
void LAPACK_EXPORT_zlaqr3(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ktop, aocl_int_t* kbot, aocl_int_t* nw, dcomplex* h, aocl_int_t* ldh, aocl_int_t* iloz, aocl_int_t* ihiz, dcomplex* z, aocl_int_t* ldz, aocl_int_t* ns, aocl_int_t* nd, dcomplex* sh, dcomplex* v, aocl_int_t* ldv, aocl_int_t* nh, dcomplex* t, aocl_int_t* ldt, aocl_int_t* nv, dcomplex* wv, aocl_int_t* ldwv, dcomplex* work, aocl_int_t* lwork);
void LAPACK_EXPORT_zlaqr4(logical* wantt, logical* wantz, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, dcomplex* h, aocl_int_t* ldh, dcomplex* w, aocl_int_t* iloz, aocl_int_t* ihiz, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zlaqr5(logical* wantt, logical* wantz, aocl_int_t* kacc22, aocl_int_t* n, aocl_int_t* ktop, aocl_int_t* kbot, aocl_int_t* nshfts, dcomplex* s, dcomplex* h, aocl_int_t* ldh, aocl_int_t* iloz, aocl_int_t* ihiz, dcomplex* z, aocl_int_t* ldz, dcomplex* v, aocl_int_t* ldv, dcomplex* u, aocl_int_t* ldu, aocl_int_t* nv, dcomplex* wv, aocl_int_t* ldwv, aocl_int_t* nh, dcomplex* wh, aocl_int_t* ldwh);
void LAPACK_EXPORT_zlaqsb(char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, double* s, double* scond, double* amax, char* equed);
void LAPACK_EXPORT_zlaqsp(char* uplo, aocl_int_t* n, dcomplex* ap, double* s, double* scond, double* amax, char* equed);
void LAPACK_EXPORT_zlaqsy(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* s, double* scond, double* amax, char* equed);
void LAPACK_EXPORT_zlar1v(aocl_int_t* n, aocl_int_t* b1, aocl_int_t* bn, double* lambda, double* d, double* l, double* ld, double* lld, double* pivmin, double* gaptol, dcomplex* z, logical* wantnc, aocl_int_t* negcnt, double* ztz, double* mingma, aocl_int_t* r, aocl_int_t* isuppz, double* nrminv, double* resid, double* rqcorr, double* work);
void LAPACK_EXPORT_zlar2v(aocl_int_t* n, dcomplex* x, dcomplex* y, dcomplex* z, aocl_int_t* incx, double* c, dcomplex* s, aocl_int_t* incc);
void LAPACK_EXPORT_zlarcm(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* c, aocl_int_t* ldc, double* rwork);
void LAPACK_EXPORT_zlarfb(char* side, char* trans, char* direct, char* storev, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* v, aocl_int_t* ldv, dcomplex* t, aocl_int_t* ldt, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* ldwork);
void LAPACK_EXPORT_zlarf(char* side, aocl_int_t* m, aocl_int_t* n, dcomplex* v, aocl_int_t* incv, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work);
void LAPACK_EXPORT_zlarfg(aocl_int_t* n, dcomplex* alpha, dcomplex* x, aocl_int_t* incx, dcomplex* tau);
void LAPACK_EXPORT_zlarfgp(aocl_int_t* n, dcomplex* alpha, dcomplex* x, aocl_int_t* incx, dcomplex* tau);
void LAPACK_EXPORT_zlarft(char* direct, char* storev, aocl_int_t* n, aocl_int_t* k, dcomplex* v, aocl_int_t* ldv, dcomplex* tau, dcomplex* t, aocl_int_t* ldt);
void LAPACK_EXPORT_zlarfx(char* side, aocl_int_t* m, aocl_int_t* n, dcomplex* v, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work);
void LAPACK_EXPORT_zlarfy(char* uplo, aocl_int_t* n, dcomplex* v, aocl_int_t* incv, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work);
void LAPACK_EXPORT_zlargv(aocl_int_t* n, dcomplex* x, aocl_int_t* incx, dcomplex* y, aocl_int_t* incy, double* c, aocl_int_t* incc);
void LAPACK_EXPORT_zlarnv(aocl_int_t* idist, aocl_int_t* iseed, aocl_int_t* n, dcomplex* x);
void LAPACK_EXPORT_zlarrv(aocl_int_t* n, double* vl, double* vu, double* d, double* l, double* pivmin, aocl_int_t* isplit, aocl_int_t* m, aocl_int_t* dol, aocl_int_t* dou, double* minrgp, double* rtol1, double* rtol2, double* w, double* werr, double* wgap, aocl_int_t* iblock, aocl_int_t* indexw, double* gers, dcomplex* z, aocl_int_t* ldz, aocl_int_t* isuppz, double* work, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_zlarscl2(aocl_int_t* m, aocl_int_t* n, double* d, dcomplex* x, aocl_int_t* ldx);
void LAPACK_EXPORT_zlartg(dcomplex* f, dcomplex* g, double* cs, dcomplex* sn, dcomplex* r__);
void LAPACK_EXPORT_zlartv(aocl_int_t* n, dcomplex* x, aocl_int_t* incx, dcomplex* y, aocl_int_t* incy, double* c, dcomplex* s, aocl_int_t* incc);
void LAPACK_EXPORT_zlarzb(char* side, char* trans, char* direct, char* storev, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, dcomplex* v, aocl_int_t* ldv, dcomplex* t, aocl_int_t* ldt, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* ldwork);
void LAPACK_EXPORT_zlarz(char* side, aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, dcomplex* v, aocl_int_t* incv, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work);
void LAPACK_EXPORT_zlarzt(char* direct, char* storev, aocl_int_t* n, aocl_int_t* k, dcomplex* v, aocl_int_t* ldv, dcomplex* tau, dcomplex* t, aocl_int_t* ldt);
void LAPACK_EXPORT_zlascl2(aocl_int_t* m, aocl_int_t* n, double* d, dcomplex* x, aocl_int_t* ldx);
void LAPACK_EXPORT_zlascl(char* type, aocl_int_t* kl, aocl_int_t* ku, double* cfrom, double* cto, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_zlaset(char* uplo, aocl_int_t* m, aocl_int_t* n, dcomplex* alpha, dcomplex* beta, dcomplex* a, aocl_int_t* lda);
void LAPACK_EXPORT_zlasr(char* side, char* pivot, char* direct, aocl_int_t* m, aocl_int_t* n, double* c, double* s, dcomplex* a, aocl_int_t* lda);
void LAPACK_EXPORT_zlassq(aocl_int_t* n, dcomplex* x, aocl_int_t* incx, double* scale, double* sumsq);
void LAPACK_EXPORT_zlaswlq(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, aocl_int_t* nb, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* ldt, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zlaswp(aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* k1, aocl_int_t* k2, aocl_int_t* ipiv, aocl_int_t* incx);
void LAPACK_EXPORT_zla_syamv(aocl_int_t* uplo, aocl_int_t* n, double* alpha, dcomplex* a, aocl_int_t* lda, dcomplex* x, aocl_int_t* incx, double* beta, double* y, aocl_int_t* incy);
void LAPACK_EXPORT_zlasyf_aa(char* uplo, aocl_int_t* j1, aocl_int_t* m, aocl_int_t* nb, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* h, aocl_int_t* ldh, dcomplex* work);
void LAPACK_EXPORT_zlasyf(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_zlasyf_rk(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, dcomplex* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_zlasyf_rook(char* uplo, aocl_int_t* n, aocl_int_t* nb, aocl_int_t* kb, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* w, aocl_int_t* ldw, aocl_int_t* info);
void LAPACK_EXPORT_zla_syrfsx_extended(aocl_int_t* prec_type, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, logical* colequ, double* c, dcomplex* b, aocl_int_t* ldb, dcomplex* y, aocl_int_t* ldy, double* berr_out, aocl_int_t* n_norms, double* err_bnds_norm, double* err_bnds_comp, dcomplex* res, double* ayb, dcomplex* dy, dcomplex* y_tail, double* rcond, aocl_int_t* ithresh, double* rthresh, double* dz_ub, logical* ignore_cwise, aocl_int_t* info);
void LAPACK_EXPORT_zlat2c(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, scomplex* sa, aocl_int_t* ldsa, aocl_int_t* info);
void LAPACK_EXPORT_zlatbs(char* uplo, char* trans, char* diag, char* normin, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, dcomplex* x, double* scale, double* cnorm, aocl_int_t* info);
void LAPACK_EXPORT_zlatdf(aocl_int_t* ijob, aocl_int_t* n, dcomplex* z, aocl_int_t* ldz, dcomplex* rhs, double* rdsum, double* rdscal, aocl_int_t* ipiv, aocl_int_t* jpiv);
void LAPACK_EXPORT_zlatps(char* uplo, char* trans, char* diag, char* normin, aocl_int_t* n, dcomplex* ap, dcomplex* x, double* scale, double* cnorm, aocl_int_t* info);
void LAPACK_EXPORT_zlatrd(char* uplo, aocl_int_t* n, aocl_int_t* nb, dcomplex* a, aocl_int_t* lda, double* e, dcomplex* tau, dcomplex* w, aocl_int_t* ldw);
void LAPACK_EXPORT_zlatrs(char* uplo, char* trans, char* diag, char* normin, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* x, double* scale, double* cnorm, aocl_int_t* info);
void LAPACK_EXPORT_zlatrz(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work);
void LAPACK_EXPORT_zlatsqr(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, aocl_int_t* nb, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* ldt, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zlatzm(char* side, aocl_int_t* m, aocl_int_t* n, dcomplex* v, aocl_int_t* incv, dcomplex* tau, dcomplex* c1, dcomplex* c2, aocl_int_t* ldc, dcomplex* work);
void LAPACK_EXPORT_zlaunhr_col_getrfnp2(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* d, aocl_int_t* info);
void LAPACK_EXPORT_zlaunhr_col_getrfnp(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* d, aocl_int_t* info);
void LAPACK_EXPORT_zlauu2(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_zlauum(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_zla_wwaddw(aocl_int_t* n, dcomplex* x, dcomplex* y, dcomplex* w);
void LAPACK_EXPORT_zpbcon(char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, double* anorm, double* rcond, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zpbequ(char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, double* s, double* scond, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_zpbrfs(char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, dcomplex* ab, aocl_int_t* ldab, dcomplex* afb, aocl_int_t* ldafb, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zpbstf(char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, aocl_int_t* info);
void LAPACK_EXPORT_zpbsv(char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, dcomplex* ab, aocl_int_t* ldab, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zpbsvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, dcomplex* ab, aocl_int_t* ldab, dcomplex* afb, aocl_int_t* ldafb, char* equed, double* s, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zpbtf2(char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, aocl_int_t* info);
void LAPACK_EXPORT_zpbtrf(char* uplo, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, aocl_int_t* info);
void LAPACK_EXPORT_zpbtrs(char* uplo, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, dcomplex* ab, aocl_int_t* ldab, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zpftrf(char* transr, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* info);
void LAPACK_EXPORT_zpftri(char* transr, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* info);
void LAPACK_EXPORT_zpftrs(char* transr, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zpocon(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* anorm, double* rcond, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zpoequb(aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* s, double* scond, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_zpoequ(aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* s, double* scond, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_zporfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zporfsx(char* uplo, char* equed, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, double* s, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zposv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zposvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, char* equed, double* s, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zposvxx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, char* equed, double* s, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* rpvgrw, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zpotf2(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_zpotrf2(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_zpotrf(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_zpotri(char* uplo, aocl_int_t* n, dcomplex* buff_A, aocl_int_t* ldim_A, aocl_int_t* info);
void LAPACK_EXPORT_zpotrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zppcon(char* uplo, aocl_int_t* n, dcomplex* ap, double* anorm, double* rcond, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zppequ(char* uplo, aocl_int_t* n, dcomplex* ap, double* s, double* scond, double* amax, aocl_int_t* info);
void LAPACK_EXPORT_zpprfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* ap, dcomplex* afp, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zppsv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* ap, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zppsvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* ap, dcomplex* afp, char* equed, double* s, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zpptrf(char* uplo, aocl_int_t* n, dcomplex* ap, aocl_int_t* info);
void LAPACK_EXPORT_zpptri(char* uplo, aocl_int_t* n, dcomplex* ap, aocl_int_t* info);
void LAPACK_EXPORT_zpptrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* ap, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zpstf2(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* piv, aocl_int_t* rank, double* tol, double* work, aocl_int_t* info);
void LAPACK_EXPORT_zpstrf(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* piv, aocl_int_t* rank, double* tol, double* work, aocl_int_t* info);
void LAPACK_EXPORT_zptcon(aocl_int_t* n, double* d, dcomplex* e, double* anorm, double* rcond, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zpteqr(char* compz, aocl_int_t* n, double* d, double* e, dcomplex* z, aocl_int_t* ldz, double* work, aocl_int_t* info);
void LAPACK_EXPORT_zptrfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* d, dcomplex* e, double* df, dcomplex* ef, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zptsv(aocl_int_t* n, aocl_int_t* nrhs, double* d, dcomplex* e, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zptsvx(char* fact, aocl_int_t* n, aocl_int_t* nrhs, double* d, dcomplex* e, double* df, dcomplex* ef, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zpttrf(aocl_int_t* n, double* d, dcomplex* e, aocl_int_t* info);
void LAPACK_EXPORT_zpttrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, double* d, dcomplex* e, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zptts2(aocl_int_t* iuplo, aocl_int_t* n, aocl_int_t* nrhs, double* d, dcomplex* e, dcomplex* b, aocl_int_t* ldb);
void LAPACK_EXPORT_zrot(aocl_int_t* n, dcomplex* cx, aocl_int_t* incx, dcomplex* cy, aocl_int_t* incy, double* c, dcomplex* s);
void LAPACK_EXPORT_zspcon(char* uplo, aocl_int_t* n, dcomplex* ap, aocl_int_t* ipiv, double* anorm, double* rcond, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zspmv(char* uplo, aocl_int_t* n, dcomplex* alpha, dcomplex* ap, dcomplex* x, aocl_int_t* incx, dcomplex* beta, dcomplex* y, aocl_int_t* incy);
void LAPACK_EXPORT_zspr(char* uplo, aocl_int_t* n, dcomplex* alpha, dcomplex* x, aocl_int_t* incx, dcomplex* ap);
void LAPACK_EXPORT_zsprfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* ap, dcomplex* afp, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zspsv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* ap, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zspsvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* ap, dcomplex* afp, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zsptrf(char* uplo, aocl_int_t* n, dcomplex* ap, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zsptri(char* uplo, aocl_int_t* n, dcomplex* ap, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zsptrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* ap, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zstedc(char* compz, aocl_int_t* n, double* d, double* e, dcomplex* z, aocl_int_t* ldz, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_zstegr(char* jobz, char* range, aocl_int_t* n, double* d, double* e, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, double* abstol, aocl_int_t* m, double* w, dcomplex* z, aocl_int_t* ldz, aocl_int_t* isuppz, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_zstein(aocl_int_t* n, double* d, double* e, aocl_int_t* m, double* w, aocl_int_t* iblock, aocl_int_t* isplit, dcomplex* z, aocl_int_t* ldz, double* work, aocl_int_t* iwork, aocl_int_t* ifail, aocl_int_t* info);
void LAPACK_EXPORT_zstemr(char* jobz, char* range, aocl_int_t* n, double* d, double* e, double* vl, double* vu, aocl_int_t* il, aocl_int_t* iu, aocl_int_t* m, double* w, dcomplex* z, aocl_int_t* ldz, aocl_int_t* nzc, aocl_int_t* isuppz, logical* tryrac, double* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_zsteqr(char* compz, aocl_int_t* n, double* d, double* e, dcomplex* z, aocl_int_t* ldz, double* work, aocl_int_t* info);
void LAPACK_EXPORT_zsteqr(char* jobz, aocl_int_t* n, double* d, double* e, dcomplex* z, aocl_int_t* ldz, double* work, aocl_int_t* info);
void LAPACK_EXPORT_zsycon_3(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, double* anorm, double* rcond, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zsycon(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, double* anorm, double* rcond, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zsycon_rook(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, double* anorm, double* rcond, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zsyconv(char* uplo, char* way, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zsyconvf(char* uplo, char* way, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zsyconvf_rook(char* uplo, char* way, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zsyequb(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* s, double* scond, double* amax, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zsymv(char* uplo, aocl_int_t* n, dcomplex* alpha, dcomplex* a, aocl_int_t* lda, dcomplex* x, aocl_int_t* incx, dcomplex* beta, dcomplex* y, aocl_int_t* incy);
void LAPACK_EXPORT_zsyr(char* uplo, aocl_int_t* n, dcomplex* alpha, dcomplex* x, aocl_int_t* incx, dcomplex* a, aocl_int_t* lda);
void LAPACK_EXPORT_zsyrfs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zsyrfsx(char* uplo, char* equed, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, double* s, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zsysv_aa_2stage(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zsysv_aa(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zsysv(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zsysv_rk(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zsysv_rook(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zsysvx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zsysvxx(char* fact, char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* af, aocl_int_t* ldaf, aocl_int_t* ipiv, char* equed, double* s, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* rcond, double* rpvgrw, double* berr, aocl_int_t* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, aocl_int_t* nparams, double* params, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_zsyswapr(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* i1, aocl_int_t* i2);
void LAPACK_EXPORT_zsytf2(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zsytf2_rk(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zsytf2_rook(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, aocl_int_t* info);
void LAPACK_EXPORT_zsytrf_aa_2stage(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zsytrf_aa(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zsytrf(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zsytrf_rk(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zsytrf_rook(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zsytri2(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zsytri2x(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* nb, aocl_int_t* info);
void LAPACK_EXPORT_zsytri_3(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zsytri_3x(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* nb, aocl_int_t* info);
void LAPACK_EXPORT_zsytri(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zsytri_rook(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zsytrs2(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zsytrs_3(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* e, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zsytrs_aa_2stage(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* tb, aocl_int_t* ltb, aocl_int_t* ipiv, aocl_int_t* ipiv2, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zsytrs_aa(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zsytrs(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_zsytrs_rook(char* uplo, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, aocl_int_t* ipiv, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_ztbcon(char* norm, char* uplo, char* diag, aocl_int_t* n, aocl_int_t* kd, dcomplex* ab, aocl_int_t* ldab, double* rcond, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ztbrfs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, dcomplex* ab, aocl_int_t* ldab, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ztbtrs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* kd, aocl_int_t* nrhs, dcomplex* ab, aocl_int_t* ldab, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_ztfsm(char* transr, char* side, char* uplo, char* trans, char* diag, aocl_int_t* m, aocl_int_t* n, dcomplex* alpha, dcomplex* a, dcomplex* b, aocl_int_t* ldb);
void LAPACK_EXPORT_ztftri(char* transr, char* uplo, char* diag, aocl_int_t* n, dcomplex* a, aocl_int_t* info);
void LAPACK_EXPORT_ztfttp(char* transr, char* uplo, aocl_int_t* n, dcomplex* arf, dcomplex* ap, aocl_int_t* info);
void LAPACK_EXPORT_ztfttr(char* transr, char* uplo, aocl_int_t* n, dcomplex* arf, dcomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_ztgevc(char* side, char* howmny, logical* select, aocl_int_t* n, dcomplex* s, aocl_int_t* lds, dcomplex* p, aocl_int_t* ldp, dcomplex* vl, aocl_int_t* ldvl, dcomplex* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ztgex2(logical* wantq, logical* wantz, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* q, aocl_int_t* ldq, dcomplex* z, aocl_int_t* ldz, aocl_int_t* j1, aocl_int_t* info);
void LAPACK_EXPORT_ztgexc(logical* wantq, logical* wantz, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* q, aocl_int_t* ldq, dcomplex* z, aocl_int_t* ldz, aocl_int_t* ifst, aocl_int_t* ilst, aocl_int_t* info);
void LAPACK_EXPORT_ztgsen(aocl_int_t* ijob, logical* wantq, logical* wantz, logical* select, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* q, aocl_int_t* ldq, dcomplex* z, aocl_int_t* ldz, aocl_int_t* m, double* pl, double* pr, double* dif, dcomplex* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* liwork, aocl_int_t* info);
void LAPACK_EXPORT_ztgsja(char* jobu, char* jobv, char* jobq, aocl_int_t* m, aocl_int_t* p, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, double* tola, double* tolb, double* alpha, double* beta, dcomplex* u, aocl_int_t* ldu, dcomplex* v, aocl_int_t* ldv, dcomplex* q, aocl_int_t* ldq, dcomplex* work, aocl_int_t* ncycle, aocl_int_t* info);
void LAPACK_EXPORT_ztgsna(char* job, char* howmny, logical* select, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* vl, aocl_int_t* ldvl, dcomplex* vr, aocl_int_t* ldvr, double* s, double* dif, aocl_int_t* mm, aocl_int_t* m, dcomplex* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_ztgsy2(char* trans, aocl_int_t* ijob, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* c, aocl_int_t* ldc, dcomplex* d, aocl_int_t* ldd, dcomplex* e, aocl_int_t* lde, dcomplex* f, aocl_int_t* ldf, double* scale, double* rdsum, double* rdscal, aocl_int_t* info);
void LAPACK_EXPORT_ztgsyl(char* trans, aocl_int_t* ijob, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* c, aocl_int_t* ldc, dcomplex* d, aocl_int_t* ldd, dcomplex* e, aocl_int_t* lde, dcomplex* f, aocl_int_t* ldf, double* scale, double* dif, dcomplex* work, aocl_int_t* lwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_ztpcon(char* norm, char* uplo, char* diag, aocl_int_t* n, dcomplex* ap, double* rcond, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ztplqt2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_ztplqt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, aocl_int_t* mb, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* t, aocl_int_t* ldt, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_ztpmlqt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, aocl_int_t* mb, dcomplex* v, aocl_int_t* ldv, dcomplex* t, aocl_int_t* ldt, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_ztpmqrt(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, aocl_int_t* nb, dcomplex* v, aocl_int_t* ldv, dcomplex* t, aocl_int_t* ldt, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_ztpqrt2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* t, aocl_int_t* ldt, aocl_int_t* info);
void LAPACK_EXPORT_ztpqrt(aocl_int_t* m, aocl_int_t* n, aocl_int_t* l, aocl_int_t* nb, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* t, aocl_int_t* ldt, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_ztprfb(char* side, char* trans, char* direct, char* storev, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, dcomplex* v, aocl_int_t* ldv, dcomplex* t, aocl_int_t* ldt, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* work, aocl_int_t* ldwork);
void LAPACK_EXPORT_ztprfs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* ap, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ztptri(char* uplo, char* diag, aocl_int_t* n, dcomplex* ap, aocl_int_t* info);
void LAPACK_EXPORT_ztptrs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* ap, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_ztpttf(char* transr, char* uplo, aocl_int_t* n, dcomplex* ap, dcomplex* arf, aocl_int_t* info);
void LAPACK_EXPORT_ztpttr(char* uplo, aocl_int_t* n, dcomplex* ap, dcomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_ztrcon(char* norm, char* uplo, char* diag, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, double* rcond, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ztrevc3(char* side, char* howmny, logical* select, aocl_int_t* n, dcomplex* t, aocl_int_t* ldt, dcomplex* vl, aocl_int_t* ldvl, dcomplex* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* info);
void LAPACK_EXPORT_ztrevc(char* side, char* howmny, logical* select, aocl_int_t* n, dcomplex* t, aocl_int_t* ldt, dcomplex* vl, aocl_int_t* ldvl, dcomplex* vr, aocl_int_t* ldvr, aocl_int_t* mm, aocl_int_t* m, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ztrexc(char* compq, aocl_int_t* n, dcomplex* t, aocl_int_t* ldt, dcomplex* q, aocl_int_t* ldq, aocl_int_t* ifst, aocl_int_t* ilst, aocl_int_t* info);
void LAPACK_EXPORT_ztrrfs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* x, aocl_int_t* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ztrsen(char* job, char* compq, logical* select, aocl_int_t* n, dcomplex* t, aocl_int_t* ldt, dcomplex* q, aocl_int_t* ldq, dcomplex* w, aocl_int_t* m, double* s, double* sep, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_ztrsna(char* job, char* howmny, logical* select, aocl_int_t* n, dcomplex* t, aocl_int_t* ldt, dcomplex* vl, aocl_int_t* ldvl, dcomplex* vr, aocl_int_t* ldvr, double* s, double* sep, aocl_int_t* mm, aocl_int_t* m, dcomplex* work, aocl_int_t* ldwork, double* rwork, aocl_int_t* info);
void LAPACK_EXPORT_ztrsyl(char* trana, char* tranb, aocl_int_t* isgn, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, dcomplex* c, aocl_int_t* ldc, double* scale, aocl_int_t* info);
void LAPACK_EXPORT_ztrti2(char* uplo, char* diag, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_ztrtri(char* uplo, char* diag, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_ztrtrs(char* uplo, char* trans, char* diag, aocl_int_t* n, aocl_int_t* nrhs, dcomplex* a, aocl_int_t* lda, dcomplex* b, aocl_int_t* ldb, aocl_int_t* info);
void LAPACK_EXPORT_ztrttf(char* transr, char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* arf, aocl_int_t* info);
void LAPACK_EXPORT_ztrttp(char* uplo, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* ap, aocl_int_t* info);
void LAPACK_EXPORT_ztzrqf(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tau, aocl_int_t* info);
void LAPACK_EXPORT_ztzrzf(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunbdb1(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, dcomplex* x11, aocl_int_t* ldx11, dcomplex* x21, aocl_int_t* ldx21, double* theta, double* phi, dcomplex* taup1, dcomplex* taup2, dcomplex* tauq1, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunbdb2(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, dcomplex* x11, aocl_int_t* ldx11, dcomplex* x21, aocl_int_t* ldx21, double* theta, double* phi, dcomplex* taup1, dcomplex* taup2, dcomplex* tauq1, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunbdb3(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, dcomplex* x11, aocl_int_t* ldx11, dcomplex* x21, aocl_int_t* ldx21, double* theta, double* phi, dcomplex* taup1, dcomplex* taup2, dcomplex* tauq1, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunbdb4(aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, dcomplex* x11, aocl_int_t* ldx11, dcomplex* x21, aocl_int_t* ldx21, double* theta, double* phi, dcomplex* taup1, dcomplex* taup2, dcomplex* tauq1, dcomplex* phantom, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunbdb5(aocl_int_t* m1, aocl_int_t* m2, aocl_int_t* n, dcomplex* x1, aocl_int_t* incx1, dcomplex* x2, aocl_int_t* incx2, dcomplex* q1, aocl_int_t* ldq1, dcomplex* q2, aocl_int_t* ldq2, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunbdb6(aocl_int_t* m1, aocl_int_t* m2, aocl_int_t* n, dcomplex* x1, aocl_int_t* incx1, dcomplex* x2, aocl_int_t* incx2, dcomplex* q1, aocl_int_t* ldq1, dcomplex* q2, aocl_int_t* ldq2, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunbdb(char* trans, char* signs, aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, dcomplex* x11, aocl_int_t* ldx11, dcomplex* x12, aocl_int_t* ldx12, dcomplex* x21, aocl_int_t* ldx21, dcomplex* x22, aocl_int_t* ldx22, double* theta, double* phi, dcomplex* taup1, dcomplex* taup2, dcomplex* tauq1, dcomplex* tauq2, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zuncsd2by1(char* jobu1, char* jobu2, char* jobv1t, aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, dcomplex* x11, aocl_int_t* ldx11, dcomplex* x21, aocl_int_t* ldx21, double* theta, dcomplex* u1, aocl_int_t* ldu1, dcomplex* u2, aocl_int_t* ldu2, dcomplex* v1t, aocl_int_t* ldv1t, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_zuncsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, char* signs, aocl_int_t* m, aocl_int_t* p, aocl_int_t* q, dcomplex* x11, aocl_int_t* ldx11, dcomplex* x12, aocl_int_t* ldx12, dcomplex* x21, aocl_int_t* ldx21, dcomplex* x22, aocl_int_t* ldx22, double* theta, dcomplex* u1, aocl_int_t* ldu1, dcomplex* u2, aocl_int_t* ldu2, dcomplex* v1t, aocl_int_t* ldv1t, dcomplex* v2t, aocl_int_t* ldv2t, dcomplex* work, aocl_int_t* lwork, double* rwork, aocl_int_t* lrwork, aocl_int_t* iwork, aocl_int_t* info);
void LAPACK_EXPORT_zung2l(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zung2r(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zungbr(char* vect, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunghr(aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zungl2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zunglq(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zungql(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zungqr(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zungr2(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zungrq(aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zungtr(char* uplo, aocl_int_t* m, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zungtsqr(aocl_int_t* m, aocl_int_t* n, aocl_int_t* mb, aocl_int_t* nb, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* ldt, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunhr_col(aocl_int_t* m, aocl_int_t* n, aocl_int_t* nb, dcomplex* a, aocl_int_t* lda, dcomplex* t, aocl_int_t* ldt, dcomplex* d, aocl_int_t* info);
void LAPACK_EXPORT_zunm22(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* n1, aocl_int_t* n2, dcomplex* q, aocl_int_t* ldq, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunm2l(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zunm2r(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zunmbr(char* vect, char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunmhr(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* ilo, aocl_int_t* ihi, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunml2(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zunmlq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunmql(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunmqr(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunmr2(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zunmr3(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zunmrq(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunmrz(char* side, char* trans, aocl_int_t* m, aocl_int_t* n, aocl_int_t* k, aocl_int_t* l, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zunmtr(char* side, char* uplo, char* trans, aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* lwork, aocl_int_t* info);
void LAPACK_EXPORT_zupgtr(char* uplo, aocl_int_t* n, dcomplex* ap, dcomplex* tau, dcomplex* q, aocl_int_t* ldq, dcomplex* work, aocl_int_t* info);
void LAPACK_EXPORT_zupmtr(char* side, char* uplo, char* trans, aocl_int_t* m, aocl_int_t* n, dcomplex* ap, dcomplex* tau, dcomplex* c, aocl_int_t* ldc, dcomplex* work, aocl_int_t* info);
logical LAPACK_EXPORT_disnan(double* din);
logical LAPACK_EXPORT_dlaisnan(double* din1, double* din2);
logical LAPACK_EXPORT_sisnan(float* sin__);
logical LAPACK_EXPORT_slaisnan(float* sin1, float* sin2);
void LAPACK_EXPORT_chla_transtype(char* ret_val, aocl_int_t* trans);
#ifdef FLA_ENABLE_VOID_RETURN_COMPLEX_FUNCTION
void LAPACK_EXPORT_cladiv(scomplex* ret_val, scomplex* x, scomplex* y);
void LAPACK_EXPORT_zladiv(dcomplex* ret_val, dcomplex* x, dcomplex* y);
#else
scomplex LAPACK_EXPORT_cladiv(scomplex* x, scomplex* y);
dcomplex LAPACK_EXPORT_zladiv(dcomplex* x, dcomplex* y);
#endif
void LAPACK_EXPORT_sgetrfnp(aocl_int_t* m, aocl_int_t* n, float* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_dgetrfnp(aocl_int_t* m, aocl_int_t* n, double* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_cgetrfnp(aocl_int_t* m, aocl_int_t* n, scomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_zgetrfnp(aocl_int_t* m, aocl_int_t* n, dcomplex* a, aocl_int_t* lda, aocl_int_t* info);
void LAPACK_EXPORT_sspffrt2( float  *ap, aocl_int_t *n, aocl_int_t * ncolm, float  *work, float  *work2 );
void LAPACK_EXPORT_dspffrt2( double *ap, aocl_int_t *n, aocl_int_t * ncolm, double *work, double *work2 );
void LAPACK_EXPORT_cspffrt2( scomplex *ap, aocl_int_t *n, aocl_int_t * ncolm, scomplex *work, scomplex *work2 );
void LAPACK_EXPORT_zspffrt2( dcomplex *ap, aocl_int_t *n, aocl_int_t * ncolm, dcomplex *work, dcomplex *work2 );
void LAPACK_EXPORT_sspffrtx( float  *ap, aocl_int_t *n, aocl_int_t * ncolm, float  *work, float  *work2 );
void LAPACK_EXPORT_dspffrtx( double *ap, aocl_int_t *n, aocl_int_t * ncolm, double *work, double *work2 );
void LAPACK_EXPORT_cspffrtx( scomplex *ap, aocl_int_t *n, aocl_int_t * ncolm, scomplex *work, scomplex *work2 );
void LAPACK_EXPORT_zspffrtx( dcomplex *ap, aocl_int_t *n, aocl_int_t * ncolm, dcomplex *work, dcomplex *work2 );
void LAPACK_EXPORT_cgetrfnpi(aocl_int_t *m, aocl_int_t *n, aocl_int_t *nfact, scomplex *a, aocl_int_t *lda, aocl_int_t *info);
void LAPACK_EXPORT_dgetrfnpi(aocl_int_t *m, aocl_int_t *n, aocl_int_t *nfact, double *a, aocl_int_t * lda, aocl_int_t *info);
void LAPACK_EXPORT_sgetrfnpi(aocl_int_t *m, aocl_int_t *n, aocl_int_t *nfact, float *a, aocl_int_t *lda, aocl_int_t *info);
void LAPACK_EXPORT_zgetrfnpi(aocl_int_t *m, aocl_int_t *n, aocl_int_t*nfact, dcomplex *a, aocl_int_t *lda, aocl_int_t *info);

void LAPACK_EXPORT_sopgtr(char *uplo, aocl_int_t *n, float *ap, float *tau, float *q, aocl_int_t *ldq, float *work, aocl_int_t *info);
void LAPACK_EXPORT_dopgtr(char *uplo, aocl_int_t *n, double *ap, double *tau, double *q, aocl_int_t *ldq, double *work, aocl_int_t *info);
void LAPACK_EXPORT_sorcsd(char *jobu1, char *jobu2, char *jobv1t, char * jobv2t, char *trans, char *signs, aocl_int_t *m, aocl_int_t *p, aocl_int_t *q, float *x11, aocl_int_t *ldx11, float *x12, aocl_int_t *ldx12, float *x21, aocl_int_t *ldx21, float *x22, aocl_int_t *ldx22, float *theta, float *u1, aocl_int_t *ldu1, float *u2, aocl_int_t *ldu2, float *v1t, aocl_int_t *ldv1t, float *v2t, aocl_int_t *ldv2t, float *work, aocl_int_t *lwork, aocl_int_t *iwork, aocl_int_t *info);
void LAPACK_EXPORT_dorcsd(char *jobu1, char *jobu2, char *jobv1t, char * jobv2t, char *trans, char *signs, aocl_int_t *m, aocl_int_t *p, aocl_int_t *q, double *x11, aocl_int_t *ldx11, double *x12, aocl_int_t *ldx12, double *x21, aocl_int_t *ldx21, double *x22, aocl_int_t *ldx22, double *theta, double *u1, aocl_int_t *ldu1, double *u2, aocl_int_t *ldu2, double *v1t, aocl_int_t *ldv1t, double *v2t, aocl_int_t *ldv2t, double *work, aocl_int_t *lwork, aocl_int_t *iwork, aocl_int_t *info);
void LAPACK_EXPORT_sorcsd2by1(char *jobu1, char *jobu2, char *jobv1t, aocl_int_t *m, aocl_int_t *p, aocl_int_t *q, float *x11, aocl_int_t *ldx11, float * x21, aocl_int_t *ldx21, float *theta, float *u1, aocl_int_t *ldu1, float *u2, aocl_int_t *ldu2, float *v1t, aocl_int_t *ldv1t, float *work, aocl_int_t *lwork, aocl_int_t *iwork, aocl_int_t *info);
void LAPACK_EXPORT_dorcsd2by1(char *jobu1, char *jobu2, char *jobv1t, aocl_int_t *m, aocl_int_t *p, aocl_int_t *q, double *x11, aocl_int_t *ldx11, double *x21, aocl_int_t *ldx21, double *theta, double *u1, aocl_int_t *ldu1, double *u2, aocl_int_t *ldu2, double *v1t, aocl_int_t *ldv1t, double *work, aocl_int_t *lwork, aocl_int_t *iwork, aocl_int_t *info);
void LAPACK_EXPORT_sorghr(aocl_int_t *n, aocl_int_t *ilo, aocl_int_t *ihi, float *a, aocl_int_t *lda, float *tau, float *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_dorghr(aocl_int_t *n, aocl_int_t *ilo, aocl_int_t *ihi, double *a, aocl_int_t *lda, double *tau, double *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_sormhr(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *ilo, aocl_int_t *ihi, float *a, aocl_int_t *lda, float *tau, float * c__, aocl_int_t *ldc, float *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_dormhr(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *ilo, aocl_int_t *ihi, double *a, aocl_int_t *lda, double * tau, double *c__, aocl_int_t *ldc, double *work, aocl_int_t *lwork, aocl_int_t *info);
void LAPACK_EXPORT_sgedmd(char *jobs, char *jobz, char *jobr, char *jobf, aocl_int_t *whtsvd, aocl_int_t *m, aocl_int_t *n, real *x, aocl_int_t *ldx, real *y, aocl_int_t *ldy, aocl_int_t *nrnk, real *tol, aocl_int_t *k, real *reig, real *imeig, real *z__, aocl_int_t *ldz, real *res, real *b, aocl_int_t *ldb, real *w, aocl_int_t *ldw, real *s, aocl_int_t *lds, real *work, aocl_int_t *lwork, aocl_int_t *iwork, aocl_int_t *liwork, aocl_int_t *info);
void LAPACK_EXPORT_dgedmd(char *jobs, char *jobz, char *jobr, char *jobf, aocl_int_t *whtsvd, aocl_int_t *m, aocl_int_t *n, doublereal *x, aocl_int_t *ldx, doublereal *y, aocl_int_t *ldy, aocl_int_t *nrnk, doublereal *tol, aocl_int_t *k, doublereal *reig, doublereal *imeig, doublereal *z__, aocl_int_t *ldz, doublereal *res, doublereal *b, aocl_int_t *ldb, doublereal *w, aocl_int_t *ldw, doublereal *s, aocl_int_t *lds, doublereal *work, aocl_int_t *lwork, aocl_int_t *iwork, aocl_int_t *liwork, aocl_int_t *info);
void LAPACK_EXPORT_cgedmd(char *jobs, char *jobz, char *jobr, char *jobf, aocl_int_t *whtsvd, aocl_int_t *m, aocl_int_t *n, scomplex *x, aocl_int_t *ldx, scomplex *y, aocl_int_t *ldy, aocl_int_t *nrnk, real *tol, aocl_int_t *k, scomplex *eigs, scomplex *z__, aocl_int_t *ldz, real *res, scomplex *b, aocl_int_t *ldb, scomplex *w, aocl_int_t *ldw, scomplex *s, aocl_int_t *lds, scomplex *zwork, aocl_int_t *lzwork, real *rwork, aocl_int_t *lrwork, aocl_int_t *iwork, aocl_int_t *liwork, aocl_int_t *info);
void LAPACK_EXPORT_zgedmd(char *jobs, char *jobz, char *jobr, char *jobf, aocl_int_t *whtsvd, aocl_int_t *m, aocl_int_t *n, dcomplex *x, aocl_int_t *ldx, dcomplex *y, aocl_int_t *ldy, aocl_int_t *nrnk, doublereal *tol, aocl_int_t *k, dcomplex *eigs, dcomplex *z__, aocl_int_t *ldz, doublereal *res, dcomplex *b, aocl_int_t *ldb, dcomplex *w, aocl_int_t *ldw, dcomplex *s, aocl_int_t *lds, dcomplex *zwork, aocl_int_t *lzwork, doublereal *rwork, aocl_int_t *lrwork, aocl_int_t *iwork, aocl_int_t *liwork, aocl_int_t *info);
void LAPACK_EXPORT_sgedmdq(char *jobs, char *jobz, char *jobr, char *jobq, char *jobt, char *jobf, aocl_int_t *whtsvd, aocl_int_t *m, aocl_int_t *n, real *f, aocl_int_t *ldf, real *x, aocl_int_t *ldx, real *y, aocl_int_t *ldy, aocl_int_t *nrnk, real *tol, aocl_int_t *k, real *reig, real *imeig, real *z__, aocl_int_t *ldz, real *res, real *b, aocl_int_t *ldb, real *v, aocl_int_t *ldv, real *s, aocl_int_t *lds, real *work, aocl_int_t *lwork, aocl_int_t *iwork, aocl_int_t *liwork, aocl_int_t *info);
void LAPACK_EXPORT_dgedmdq(char *jobs, char *jobz, char *jobr, char *jobq, char *jobt, char *jobf, aocl_int_t *whtsvd, aocl_int_t *m, aocl_int_t *n, doublereal *f, aocl_int_t *ldf, doublereal *x, aocl_int_t *ldx, doublereal *y, aocl_int_t *ldy, aocl_int_t *nrnk, doublereal *tol, aocl_int_t *k, doublereal *reig, doublereal *imeig, doublereal *z__, aocl_int_t *ldz, doublereal *res, doublereal *b, aocl_int_t *ldb, doublereal *v, aocl_int_t *ldv, doublereal *s, aocl_int_t *lds, doublereal *work, aocl_int_t *lwork, aocl_int_t *iwork, aocl_int_t *liwork, aocl_int_t *info);
void LAPACK_EXPORT_cgedmdq(char *jobs, char *jobz, char *jobr, char *jobq, char *jobt, char *jobf, aocl_int_t *whtsvd, aocl_int_t *m, aocl_int_t *n, scomplex *f, aocl_int_t *ldf, scomplex *x, aocl_int_t *ldx, scomplex *y, aocl_int_t *ldy, aocl_int_t *nrnk, real *tol, aocl_int_t *k, scomplex *eigs, scomplex *z__, aocl_int_t *ldz, real *res, scomplex *b, aocl_int_t *ldb, scomplex *v, aocl_int_t *ldv, scomplex *s, aocl_int_t *lds, scomplex *zwork, aocl_int_t *lzwork, real *work, aocl_int_t *lwork, aocl_int_t *iwork, aocl_int_t *liwork, aocl_int_t *info);
void LAPACK_EXPORT_zgedmdq(char *jobs, char *jobz, char *jobr, char *jobq, char *jobt, char *jobf, aocl_int_t *whtsvd, aocl_int_t *m, aocl_int_t *n, dcomplex *f, aocl_int_t *ldf, dcomplex *x, aocl_int_t *ldx, dcomplex *y, aocl_int_t *ldy, aocl_int_t *nrnk, doublereal *tol, aocl_int_t *k, dcomplex *eigs, dcomplex *z__, aocl_int_t *ldz, doublereal *res, dcomplex *b, aocl_int_t *ldb, dcomplex *v, aocl_int_t *ldv, dcomplex *s, aocl_int_t *lds, dcomplex *zwork, aocl_int_t *lzwork, doublereal *work, aocl_int_t *lwork, aocl_int_t *iwork, aocl_int_t *liwork, aocl_int_t *info);