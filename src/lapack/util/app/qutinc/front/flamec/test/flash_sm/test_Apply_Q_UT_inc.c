/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


#define N_PARAM_COMBOS    1

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

char* pc_str[N_PARAM_COMBOS] = { "" };

void time_Apply_Q_UT_inc(
               integer param_combo, integer type, integer nrepeats, integer m, integer n,
               FLA_Obj A, FLA_Obj TW, FLA_Obj W1, FLA_Obj B, FLA_Obj B_ref,
               FLA_Obj A_flat, FLA_Obj T_flat, FLA_Obj W_flat, FLA_Obj B_flat,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  integer 
    datatype,
    n_threads,
    m_input, n_input,
    m, n,
    p_first, p_last, p_inc,
    p,
    n_repeats,
    param_combo,
    i,
    n_param_combos = N_PARAM_COMBOS;
  
  fla_dim_t
    nb_flash,
    nb_alg;

  char *colors = "brkgmcbrkgmcbrkgmc";
  char *ticks  = "o+*xso+*xso+*xso+*xs";
  char m_dim_desc[14];
  char m_dim_tag[10];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj A_flat, T_flat, W_flat, B_flat;
  FLA_Obj A, TW, W1, B, B_ref;
  
  FLA_Init( );


  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter algorithmic blocksize: ", '%' );
  scanf( "%u", &nb_alg );
  fprintf( stdout, "%c %u\n", '%', nb_alg );

  fprintf( stdout, "%c enter FLASH blocksize: ", '%' );
  scanf( "%u", &nb_flash );
  fprintf( stdout, "%c %u\n", '%', nb_flash );

  fprintf( stdout, "%c enter problem size first, last, inc: ", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m n (-1 means bind to problem size): ", '%' );
  scanf( "%d%d", &m_input, &n_input );
  fprintf( stdout, "%c %d %d\n", '%', m_input, n_input );

  fprintf( stdout, "%c enter the number of SuperMatrix threads: ", '%' );
  scanf( "%d", &n_threads );
  fprintf( stdout, "%c %d\n", '%', n_threads );


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

  //datatype = FLA_FLOAT;
  //datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  datatype = FLA_DOUBLE_COMPLEX;

  FLASH_Queue_set_num_threads( n_threads );

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;
    n = n_input;

    if( m < 0 ) m = p / f2c_abs(m_input);
    if( n < 0 ) n = p / f2c_abs(n_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ )
    {
      FLA_Obj_create( datatype, m,      m, &A_flat );
      FLA_Obj_create( datatype, nb_alg, m, &T_flat );
      FLA_Obj_create( datatype, nb_alg, n, &W_flat );
      FLA_Obj_create( datatype, m,      n, &B_flat );
      FLA_Random_matrix( A_flat );
      FLA_Random_matrix( B_flat );

      FLASH_Obj_create_ext( datatype, m, n, 1, &nb_flash, &nb_flash, &B );
      FLASH_Obj_create_ext( datatype, m, n, 1, &nb_flash, &nb_flash, &B_ref );

      FLASH_QR_UT_inc_create_hier_matrices( A_flat, 1, &nb_flash, nb_alg, &A, &TW );
      FLASH_Apply_Q_UT_inc_create_workspace( TW, B, &W1 );

      FLASH_Obj_hierarchify( B_flat, B );


      fprintf( stdout, "data_qrutinc_%s( %d, 1:5 ) = [ %d  ", pc_str[param_combo], i, p );
      fflush( stdout );

      time_Apply_Q_UT_inc( param_combo, FLA_ALG_REFERENCE, n_repeats, m, n,
                           A, TW, W1, B, B_ref, A_flat, T_flat, W_flat, B_flat,
                           &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Apply_Q_UT_inc( param_combo, FLA_ALG_FRONT, n_repeats, m, n,
                           A, TW, W1, B, B_ref, A_flat, T_flat, W_flat, B_flat,
                           &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );


      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &A_flat );
      FLA_Obj_free( &T_flat );
      FLA_Obj_free( &W_flat );
      FLA_Obj_free( &B_flat );

      FLASH_Obj_free( &A );
      FLASH_Obj_free( &TW );
      FLASH_Obj_free( &W1 );
      FLASH_Obj_free( &B );
      FLASH_Obj_free( &B_ref );
    }

    fprintf( stdout, "\n" );
  }

  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_param_combos; i++ ) {
    fprintf( stdout, "plot( data_qrutinc_%s( :,1 ), data_qrutinc_%s( :, 2 ), '%c:%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_qrutinc_%s( :,1 ), data_qrutinc_%s( :, 4 ), '%c-.%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_param_combos; i++ )
    fprintf( stdout, "'ref\\_qrutinc\\_%s', 'fla\\_qrutinc\\_%s', ... \n", pc_str[i], pc_str[i] );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );


  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME qrutinc front-end performance (%s)' );\n",
           m_dim_desc );
  fprintf( stdout, "print -depsc qrutinc_front_%s.eps\n", m_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );

  FLA_Finalize( );

  return 0;
}

