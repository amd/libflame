%
%
%   Copyright (C) 2014, The University of Texas at Austin
%
%   This file is part of libflame and is available under the 3-Clause
%   BSD license, which can be found in the LICENSE file at the top-level
%   directory, or at http://opensource.org/licenses/BSD-3-Clause
%
%

function [ C_out ] = FLA_Symm_ru_blk_var8( A, B, C, nb_alg )

  [ ATL, ATR, ...
    ABL, ABR ] = FLA_Part_2x2( A, ...
                               0, 0, 'FLA_BR' );

  [ BL, BR ] = FLA_Part_1x2( B, ...
                               0, 'FLA_RIGHT' );

  [ CL, CR ] = FLA_Part_1x2( C, ...
                               0, 'FLA_RIGHT' );

  while ( size( ABR, 1 ) < size( A, 1 ) )

    b = fla_min( size( ATL, 1 ), nb_alg );

    [ A00, A01, A02, ...
      A10, A11, A12, ...
      A20, A21, A22 ] = FLA_Repart_2x2_to_3x3( ATL, ATR, ...
                                               ABL, ABR, ...
                                               b, b, 'FLA_TL' );

    [ B0, B1, B2 ]= FLA_Repart_1x2_to_1x3( BL, BR, ...
                                         b, 'FLA_LEFT' );

    [ C0, C1, C2 ]= FLA_Repart_1x2_to_1x3( CL, CR, ...
                                         b, 'FLA_LEFT' );

    %------------------------------------------------------------%

    C0 = C0 + B1 * A01';
    C1 = C1 + B1 * A11;
    C2 = C2 + B1 * A12;

    %------------------------------------------------------------%

    [ ATL, ATR, ...
      ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00, A01, A02, ...
                                             A10, A11, A12, ...
                                             A20, A21, A22, ...
                                             'FLA_BR' );

    [ BL, BR ] = FLA_Cont_with_1x3_to_1x2( B0, B1, B2, ...
                                           'FLA_RIGHT' );

    [ CL, CR ] = FLA_Cont_with_1x3_to_1x2( C0, C1, C2, ...
                                           'FLA_RIGHT' );

  end

  C_out = [ CL, CR ];

return

