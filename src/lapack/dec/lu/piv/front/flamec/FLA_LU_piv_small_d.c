/*
    Copyright (c) 2021-2023 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif

/*
 * LU with partial pivoting for tiny matrices
 *
 * All the computations are done inline without using
 * corresponding BLAS APIs to reduce function overheads.
 */
fla_dim_t FLA_LU_piv_small_d_var0( fla_dim_t *m, fla_dim_t *n,
                                   doublereal *a, fla_dim_t *lda,
                                   aocl_int_t *ipiv,
                                   fla_dim_t *info)
{
    fla_dim_t mi, ni;
    fla_dim_t i, j, i_1;

    doublereal p_val, max_val, t_val;
    doublereal *acur, *apiv, *asrc;
    fla_dim_t p_idx;
    fla_dim_t min_m_n = fla_min(*m, *n);

    for( i = 0; i < min_m_n; i++ )
    {
        mi = *m - i;
        ni = *n - i;

        acur = &a[i + *lda * i];

        /* Find the pivot element */
        max_val = 0;
        p_idx = i;
        for( i_1 = 0; i_1 < mi; i_1++ )
        {
            t_val = acur[i_1];
            t_val = ( t_val < 0.0 ) ? -t_val : t_val;
            if( t_val > max_val )
            {
                max_val = t_val;
                p_idx = i + i_1;
            }
        }

        apiv = a + p_idx;
        asrc = a + i;
        ipiv[i] = (aocl_int_t)(p_idx + 1);

        /* Swap rows and calculate a column of L */
        if( max_val != 0.0 )
        {
            /* Swap entire rows */
            if( p_idx != i)
            {
                for( i_1 = 0; i_1 < *n; i_1++ )
                {
                    t_val = apiv[i_1 * *lda];
                    apiv[i_1 * *lda] = asrc[i_1 * *lda];
                    asrc[i_1 * *lda] = t_val;
                }
            }

            /* Calculate scalefactors (L)  & update trailing matrix */
            p_val = *acur;
            for( i_1 = 1; i_1 < mi; i_1++ )
            {
                acur[i_1] = acur[i_1] / p_val;
                for( j = 1; j < ni; j++ )
                {
                    acur[i_1 + j * *lda] = acur[i_1 + j * *lda] - acur[j * *lda] * acur[i_1];
                }
            }
        }
        else
        {
            *info = ( *info == 0 ) ? p_idx + 1 : *info;
        }
    }
    

    return *info;
}

/*
 * LU with partial pivoting for small matrices
 *
 * This is an unblocked variant making use of BLAS APIs
 */
fla_dim_t FLA_LU_piv_small_d_var1( fla_dim_t *m, fla_dim_t *n,
                                   doublereal *a, fla_dim_t *lda,
                                   aocl_int_t *ipiv,
                                   fla_dim_t *info)
{
    fla_dim_t a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    fla_dim_t c__1 = 1;
    doublereal c_n1 = -1.;


    /* Local variables */
    fla_dim_t i__, j, jp;
    extern doublereal dlamch_(char *);
    doublereal sfmin;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;

/*     Compute machine safe minimum */
    sfmin = dlamch_("S");

    i__1 = fla_min(*m,*n);
    for ( j = 1; j <= i__1; ++j )
    {

/*        Find pivot and test for singularity. */
	    i__2 = *m - j + 1;
	    jp = j - 1 + aocl_blas_idamax(&i__2, &a[j + j * a_dim1], &c__1);
	    ipiv[j] = (aocl_int_t)jp;
	    if (a[jp + j * a_dim1] != 0.)
        {

/*           Apply the interchange to columns 1:N. */
            
            if (jp != j)
            {
                aocl_blas_dswap(n, &a[j + a_dim1], lda, &a[jp + a_dim1], lda);
            }

/*           Compute elements J+1:M of J-th column. */

	        if (j < *m)
            {
                d__1 = a[j + j * a_dim1];
                d__1 = (d__1 < 0) ? -d__1 : d__1;
                if (d__1 >= sfmin)
                {
		            i__2 = *m - j;
		            d__1 = 1. / a[j + j * a_dim1];
		            aocl_blas_dscal(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
		        }
                else
                {
		            i__2 = *m - j;
		            for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        a[j + i__ + j * a_dim1] /= a[j + j * a_dim1];
/* L20: */
		             }
                 }
            }
        }
        else if ( *info == 0 )
        {
            *info = j;
        }
        
        if ( j < fla_min( *m, *n ) )
        {

/*           Update trailing submatrix. */

	        i__2 = *m - j;
            i__3 = *n - j;
            aocl_blas_dger(&i__2, &i__3, &c_n1, &a[j + 1 + j * a_dim1], &c__1, &a[j + (
		           j + 1) * a_dim1], lda, &a[j + 1 + (j + 1) * a_dim1], lda);
        }
/* L10: */
    }
    return *info;
}

/*
 * LU with partial pivoting for medium-sized matrices
 *
 * This is a simple non-recursive blocked variant making
 * use of BLAS APIs.
 */
fla_dim_t FLA_LU_piv_small_d_var2( fla_dim_t *m, fla_dim_t *n,
                                   doublereal *a, fla_dim_t *lda,
                                   aocl_int_t *ipiv,
                                   fla_dim_t *info)
{
    fla_dim_t c__1 = 1;
    doublereal c_b16 = 1.;
    doublereal c_b19 = -1.;

    fla_dim_t a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    fla_dim_t i__, j, jb, nb;
    extern doublereal dlamch_(char *);
    fla_dim_t iinfo;

#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    nb = FLA_SMALL_LU_BLOCKSIZE;
	i__1 = fla_min(*m,*n);
	i__2 = nb;
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
    {
/* Computing MIN */
	    i__3 = fla_min(*m,*n) - j + 1;
	    jb = fla_min(i__3,nb);

/*           Factor diagonal and subdiagonal blocks and test for exact
             singularity. */

	    i__3 = *m - j + 1;
	    aocl_lapack_dgetrf2(&i__3, &jb, &a_ref(j, j), lda, &ipiv[j], &iinfo);

/*           Adjust INFO and the pivot indices. */

	    if (*info == 0 && iinfo > 0)
        {
            *info = iinfo + j - 1;
	    }
/* Computing MIN */
	    i__4 = *m, i__5 = j + jb - 1;
	    i__3 = fla_min(i__4,i__5);
	    for (i__ = j; i__ <= i__3; ++i__)
        {
            ipiv[i__] = (aocl_int_t)(j - 1 + ipiv[i__]);
/* L10: */
	    }

/*           Apply interchanges to columns 1:J-1. */

	    i__3 = j - 1;
	    i__4 = j + jb - 1;
	    aocl_lapack_dlaswp(&i__3, &a[a_offset], lda, &j, &i__4, &ipiv[1], &c__1);

	    if (j + jb <= *n)
        {

/*              Apply interchanges to columns J+JB:N. */

		    i__3 = *n - j - jb + 1;
    		i__4 = j + jb - 1;
	    	aocl_lapack_dlaswp(&i__3, &a_ref(1, j + jb), lda, &j, &i__4, &ipiv[1], &
		        	c__1);

/*              Compute block row of U. */
            
            i__3 = *n - j - jb + 1;
            aocl_blas_dtrsm("Left", "Lower", "No transpose", "Unit", &jb, &i__3, &c_b16,
                    &a_ref(j, j), lda, &a_ref(j, j + jb), lda);
            if (j + jb <= *m)
            {
/*                 Update trailing submatrix. */
                i__3 = *m - j - jb + 1;
                i__4 = *n - j - jb + 1;
                aocl_blas_dgemm("No transpose", "No transpose", &i__3, &i__4, &jb,
                        &c_b19, &a_ref(j + jb, j), lda, &a_ref(j, j + jb),
                        lda, &c_b16, &a_ref(j + jb, j + jb), lda);
            }
        }
/* L20: */
	}
#undef a_ref
	return *info;
}

