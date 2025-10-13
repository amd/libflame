/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern int fla_thread_get_num_threads(void);
extern TLS_CLASS_SPEC fla_lu_t* fla_lu_piv_cntl;
extern TLS_CLASS_SPEC fla_lu_t* fla_lu_piv_cntl2;

FLA_Error FLA_LU_piv( FLA_Obj A, FLA_Obj p )
{
  FLA_Error r_val = FLA_SUCCESS;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_LU_piv_check( A, p );

  // Invoke FLA_LU_piv_internal() with large control tree.
  r_val = FLA_LU_piv_internal( A, p, fla_lu_piv_cntl2 );

  // This is invalid as FLA_LU_piv_internal returns a null pivot index.
  // Check for singularity.
  //if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
  //  r_val = FLA_LU_find_zero_on_diagonal( A );

  return r_val;
}


void FLA_get_optimum_params_getrf(fla_dim_t m, fla_dim_t n, fla_dim_t *nb, int *n_threads)
{
    int available_n_threads = fla_thread_get_num_threads();

#ifdef FLA_OPENMP_MULTITHREADING

    double ratio = (double)m / (double)n;

    // choose n_threads = 64/32/16/<16
    *n_threads = (available_n_threads >= 64   ? 64
                  : available_n_threads >= 32 ? 32
                  : available_n_threads >= 16 ? 16
                  : available_n_threads >= 4  ? available_n_threads
                                              : 1);

    // dynamic block-size: aim for ~2*threads blocks, clamp [16,128], round to 8
    {
        int dim = (m > n ? m : n);
        int blocks = *n_threads * 2;
        int nb_dyn = ((dim + blocks - 1) / blocks + 7) / 8 * 8;
        *nb = nb_dyn < 16 ? 16 : (nb_dyn > 128 ? 128 : nb_dyn);
    }

    // now override for special shapes / sizes if desired
    if(*n_threads == 64)
    {
        if(ratio > TALL_RATIO_THRESHOLD)
        {
            if(m > 8192)
            {
                *nb = 128;
                *n_threads = 48;
            }
            else if(m > 4096)
            {
                *nb = 64;
                *n_threads = 32;
            }
            else if(m > 2048)
            {
                *nb = 32;
                *n_threads = 16;
            }
            else
            {
                *nb = 16;
                *n_threads = 8;
            }
        }
        else if(ratio < WIDE_RATIO_THRESHOLD)
        {
            *nb = fla_max(*nb, (n > 8192 ? 96 : 64));
            *n_threads = 32;
        }
        else
        {
            if(m > 8192)
            {
                *nb = 128;
                *n_threads = 64;
            }
            else if(m > 4096)
            {
                *nb = 64;
                *n_threads = 32;
            }
            else if(m > 2048)
            {
                *nb = 64;
                *n_threads = 24;
            }
            else
            {
                *nb = 16;
                *n_threads = 8;
            }
        }
    }
    else if(*n_threads == 32)
    {
        if(ratio > TALL_RATIO_THRESHOLD)
        {
            *nb = fla_max(*nb, (m > 4096 ? 64 : 32));
            *n_threads = 16;
        }
        else if(ratio < WIDE_RATIO_THRESHOLD)
        {
            *nb = fla_max(*nb, (n > 4096 ? 64 : 48));
            *n_threads = 16;
        }
        else
        {
            if(m > 4096)
            {
                *nb = 128;
                *n_threads = 32;
            }
            else if(m > 2048)
            {
                *nb = 64;
                *n_threads = 16;
            }
            else
            {
                *nb = 32;
                *n_threads = 8;
            }
        }
    }
    else if(*n_threads == 16)
    {
        if(ratio > TALL_RATIO_THRESHOLD)
        {
            *nb = fla_max(*nb, (m > 2048 ? 64 : 32));
            *n_threads = 12;
        }
        else if(ratio < WIDE_RATIO_THRESHOLD)
        {
            *nb = fla_max(*nb, (n > 2048 ? 64 : 32));
            *n_threads = 12;
        }
        else
        {
            *nb = fla_max(*nb, (m > 4096 ? 128 : (m > 1024 ? 64 : 32)));
            *n_threads = 16;
        }
    }
    else
    {
        // fallback for <16 threads
        if(m <= 512 || n <= 512)
        {
            *nb = fla_max(*nb, 16);
            *n_threads = 4;
        }
        else if(m <= 2048 || n <= 2048)
        {
            *nb = fla_max(*nb, 32);
            *n_threads = 8;
        }
        else
        {
            *nb = fla_max(*nb, 64); /* keep *n_threads */
        }
    }

    // never exceed hw threads
    if(*n_threads > available_n_threads)
        *n_threads = available_n_threads;
    // Round block-size up to a multiple of 4
    *nb = ((*nb + 3) / 4) * 4;

#else
    *nb = 64;
    *n_threads = 1;
#endif

    return;
}