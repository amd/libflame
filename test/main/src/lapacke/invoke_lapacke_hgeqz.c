/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_hgeqz(integer datatype, int layout, char job, char compq, char compz,
                             integer n, integer ilo, integer ihi, void *h, integer ldh, void *t,
                             integer ldt, void *alpha, void *alphar, void *alphai, void *beta,
                             void *q, integer ldq, void *z, integer ldz)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_shgeqz(layout, job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar,
                                  alphai, beta, q, ldq, z, ldz);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dhgeqz(layout, job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar,
                                  alphai, beta, q, ldq, z, ldz);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_chgeqz(layout, job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha,
                                  beta, q, ldq, z, ldz);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zhgeqz(layout, job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha,
                                  beta, q, ldq, z, ldz);
            break;
        }
    }
    return info;
}
