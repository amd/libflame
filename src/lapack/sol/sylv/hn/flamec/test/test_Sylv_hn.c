/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_BLOCKED   1
#define FLA_ALG_UNBLOCKED 2
#define FLA_ALG_UNB_OPT   3


void time_Sylv_hn(
               integer variant, integer type, integer n_repeats, integer m, integer n, integer nb_alg,
               FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref, FLA_Obj scale,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  integer 
    m_input, n_input,
    m, n,
    p_first, p_last, p_inc,
    p,
    nb_alg,
    variant,
    n_repeats,
    i, j,
    datatype,
    n_variants = 18;

  integer  sign;
  
  integer  blocksize[16];

  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char n_dim_desc[14];
  char m_dim_tag[10];
  char n_dim_tag[10];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    A, B, C, C_ref, scale, isgn, norm;
  

  /* Initialize FLAME */
  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c Enter blocking size:", '%' );
  scanf( "%d", &nb_alg );
  fprintf( stdout, "%c %d\n", '%', nb_alg );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c Enter sign (-1 or 1):", '%' );
  scanf( "%d", &sign );
  fprintf( stdout, "%c %d\n", '%', sign );

  fprintf( stdout, "%c enter m n (-1 means bind to problem size): ", '%' );
  scanf( "%d %d", &m_input, &n_input );
  fprintf( stdout, "%c %d %d\n", '%', m_input, n_input );


  /* Delete all existing data structures */
  fprintf( stdout, "\nclear all;\n\n" );


  if     ( m_input >  0 ) {
    sprintf( m_dim_desc, "m = %d", m_input );
    sprintf( m_dim_tag,  "m%dc", m_input);
  }
  else if( m_input <  -1 ) {
    sprintf( m_dim_desc, "m = p/%d", -m_input );
    sprintf( m_dim_tag,  "m%dp", -m_input );
  }
  else if( m_input == -1 ) {
    sprintf( m_dim_desc, "m = p" );
    sprintf( m_dim_tag,  "m%dp", 1 );
  }
  if     ( n_input >  0 ) {
    sprintf( n_dim_desc, "n = %d", n_input );
    sprintf( n_dim_tag,  "n%dc", n_input);
  }
  else if( n_input <  -1 ) {
    sprintf( n_dim_desc, "n = p/%d", -n_input );
    sprintf( n_dim_tag,  "n%dp", -n_input );
  }
  else if( n_input == -1 ) {
    sprintf( n_dim_desc, "n = p" );
    sprintf( n_dim_tag,  "n%dp", 1 );
  }

  if ( 0 < sign )
    isgn = FLA_ONE;
  else
    isgn = FLA_MINUS_ONE;


  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {

    m = m_input;
    n = n_input;

    if( m < 0 ) m = p / f2c_abs(m_input);
    if( n < 0 ) n = p / f2c_abs(n_input);

    //datatype = FLA_FLOAT;
    //datatype = FLA_DOUBLE;
    //datatype = FLA_COMPLEX;
    datatype = FLA_DOUBLE_COMPLEX;

    FLA_Obj_create( datatype, m, m, &A );
    FLA_Obj_create( datatype, n, n, &B );
    FLA_Obj_create( datatype, m, n, &C );
    FLA_Obj_create( datatype, m, n, &C_ref );

    if ( datatype == FLA_DOUBLE || datatype == FLA_DOUBLE_COMPLEX )
    {
      FLA_Obj_create( FLA_DOUBLE, 1, 1, &scale );
      FLA_Obj_create( FLA_DOUBLE, 1, 1, &norm );
    }
    else if ( datatype == FLA_FLOAT || datatype == FLA_COMPLEX )
    {
      FLA_Obj_create( FLA_FLOAT, 1, 1, &scale );
      FLA_Obj_create( FLA_FLOAT, 1, 1, &norm );
    }

    FLA_Random_tri_matrix( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
    FLA_Random_tri_matrix( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, B );
    FLA_Random_matrix( C );

    FLA_Norm1( A, norm );
    FLA_Shift_diag( FLA_NO_CONJUGATE, norm, A );

    FLA_Norm1( B, norm );
    if ( FLA_Obj_is( isgn, FLA_MINUS_ONE ) )
      FLA_Negate( norm );
    FLA_Shift_diag( FLA_NO_CONJUGATE, norm, B );


    time_Sylv_hn( 0, FLA_ALG_REFERENCE, n_repeats, m, n, nb_alg,
                  isgn, A, B, C, C_ref, scale, &dtime, &diff, &gflops );

    fprintf( stdout, "data_REF( %d, 1:2 ) = [ %d  %6.3lf ]; \n", i, p, gflops );
    fflush( stdout );

    for ( variant = 1; variant <= n_variants; variant++ ){
      
      fprintf( stdout, "data_var%d( %d, 1:3 ) = [ %d  ", variant, i, p );
      fflush( stdout );

      time_Sylv_hn( variant, FLA_ALG_UNB_OPT, n_repeats, m, n, nb_alg,
                    isgn, A, B, C, C_ref, scale, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Sylv_hn( variant, FLA_ALG_BLOCKED, n_repeats, m, n, nb_alg,
                    isgn, A, B, C, C_ref, scale, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, " ]; \n" );
      fflush( stdout );
    }

    FLA_Obj_free( &A );
    FLA_Obj_free( &B );
    FLA_Obj_free( &C );
    FLA_Obj_free( &C_ref );
    FLA_Obj_free( &scale );
    FLA_Obj_free( &norm );
    fprintf( stdout, "\n" );
  }
/*
  // Print the MATLAB commands to plot the data

  // Delete all existing figures
  fprintf( stdout, "figure;\n" );

  // Plot the performance of the reference implementation
  fprintf( stdout, "plot( data_REF( :,1 ), data_REF( :, 2 ), '-' ); \n" );

  // Indicate that you want to add to the existing plot
  fprintf( stdout, "hold on;\n" );

  // Plot the data for the other numbers of threads
  for ( i = 1; i <= n_variants; i++ ){
    fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 2 ), '%c:%c' ); \n", 
             i, i, colors[ i-1 ], ticks[ i-1 ] );
  }

  fprintf( stdout, "legend( ... \n" );
  fprintf( stdout, "'Reference', ... \n" );

  for ( i = 1; i <= n_variants; i++ )
    fprintf( stdout, "'FLAME var%d', ... \n", i );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME sylv\\_tn performance (%s)' );\n", 
           m_dim_desc );
  fprintf( stdout, "print -depsc sylv_tn_%s.eps\n", m_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );
}

