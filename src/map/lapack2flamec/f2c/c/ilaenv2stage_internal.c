#include "FLA_f2c.h" /* > \brief \b ILAENV2STAGE */

aocl_int64_t aocl_lapack_ilaenv2stage(aocl_int64_t *ispec, char *name__, char *opts,
                                      aocl_int64_t *n1, aocl_int64_t *n2, aocl_int64_t *n3,
                                      aocl_int64_t *n4)
{
    /* System generated locals */
    aocl_int64_t ret_val;
    /* Local variables */
    aocl_int64_t iispec;
    /* -- LAPACK auxiliary routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* July 2017 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    switch(*ispec)
    {
        case 1:
            goto L10;
        case 2:
            goto L10;
        case 3:
            goto L10;
        case 4:
            goto L10;
        case 5:
            goto L10;
    }
    /* Invalid value for ISPEC */
    ret_val = -1;
    return ret_val;
L10: /* 2stage eigenvalues and SVD or related subroutines. */
    iispec = *ispec + 16;
    ret_val = aocl_lapack_iparam2stage(&iispec, name__, opts, n1, n2, n3, n4); // , name_len, opts_len);
    return ret_val;
    /* End of ILAENV2STAGE */
}
/* ilaenv2stage_ */
