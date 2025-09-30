/* maxloc.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
    on Microsoft Windows system, link with libf2c.lib;
    on Linux or Unix systems, link with .../path/to/libf2c.a -lm
    or, if you install libf2c.a in a standard place, with -lf2c -lm
    -- in that order, at the end of the command line, as in
        cc *.o -lf2c -lm
    Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

        http://www.netlib.org/f2c/libf2c.zip
*/

#include "FLA_f2c.h"

/* ********************************************************************************** */
/** Generated wrapper function */
aocl_int_t smaxloc_(real *a, aocl_int_t *dimm)
{
#if FLA_ENABLE_ILP64
    return aocl_lapack_smaxloc(a, dimm);
#else
    aocl_int64_t dimm_64 = *dimm;

    return aocl_lapack_smaxloc(a, &dimm_64);
#endif
}

aocl_int64_t aocl_lapack_smaxloc(real *a, aocl_int64_t *dimm)
{
    /* System generated locals */
    aocl_int64_t ret_val, i__1;

    /* Local variables */
    aocl_int64_t i__;
    real smax;

    /* Parameter adjustments */
    --a;

    /* Function Body */
    ret_val = 1;
    smax = a[1];
    i__1 = *dimm;
    for(i__ = 2; i__ <= i__1; ++i__)
    {
        if(smax < a[i__])
        {
            smax = a[i__];
            ret_val = i__;
        }
        /* L10: */
    }
    return ret_val;
} /* smaxloc_ */

/* ********************************************************************************** */
/** Generated wrapper function */
aocl_int_t dmaxloc_(doublereal *a, aocl_int_t *dimm)
{
#if FLA_ENABLE_ILP64
    return aocl_lapack_dmaxloc(a, dimm);
#else
    aocl_int64_t dimm_64 = *dimm;

    return aocl_lapack_dmaxloc(a, &dimm_64);
#endif
}

aocl_int64_t aocl_lapack_dmaxloc(doublereal *a, aocl_int64_t *dimm)
{
    /* System generated locals */
    aocl_int64_t ret_val, i__1;

    /* Local variables */
    aocl_int64_t i__;
    doublereal dmax__;

    /* Parameter adjustments */
    --a;

    /* Function Body */
    ret_val = 1;
    dmax__ = a[1];
    i__1 = *dimm;
    for(i__ = 2; i__ <= i__1; ++i__)
    {
        if(dmax__ < a[i__])
        {
            dmax__ = a[i__];
            ret_val = i__;
        }
        /* L20: */
    }
    return ret_val;
} /* dmaxloc_ */
