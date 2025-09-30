/* ilaver.f -- translated by f2c (version 20061008).
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

void ilaver_(aocl_int_t *vers_major__, aocl_int_t *vers_minor__, aocl_int_t *vers_patch__)
{
#if FLA_ENABLE_ILP64
    return aocl_lapack_ilaver(vers_major__, vers_minor__, vers_patch__);
#else
    aocl_int64_t vers_major___64 = *vers_major__;
    aocl_int64_t vers_minor___64 = *vers_minor__;
    aocl_int64_t vers_patch___64 = *vers_patch__;

    aocl_lapack_ilaver(&vers_major___64, &vers_minor___64, &vers_patch___64);

    *vers_major__ = (aocl_int_t)vers_major___64;
    *vers_minor__ = (aocl_int_t)vers_minor___64;
    *vers_patch__ = (aocl_int_t)vers_patch___64;
#endif
}


/* Subroutine */ void aocl_lapack_ilaver(aocl_int64_t *vers_major__, aocl_int64_t *vers_minor__, aocl_int64_t *vers_patch__)
{

    /*  -- LAPACK routine (version 3.4.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2010 */

    /*     .. Scalar Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  This subroutine return the Lapack version. */

    /*  Arguments */
    /*  ========= */

    /*  VERS_MAJOR   (output) INTEGER */
    /*      return the lapack major version */
    /*  VERS_MINOR   (output) INTEGER */
    /*      return the lapack minor version from the major version */
    /*  VERS_PATCH   (output) INTEGER */
    /*      return the lapack patch version from the minor version */

    /*     .. Executable Statements .. */
    /* Current version is 3.9.0 */
    *vers_major__ = 3;
    *vers_minor__ = 9;
    *vers_patch__ = 0;
    /*  ===================================================================== */
} /* ilaver_ */
