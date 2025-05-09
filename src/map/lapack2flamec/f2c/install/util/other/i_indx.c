#include "f2c.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers
integer i_indx(a, b, la, lb) char *a, *b;
ftnlen la, lb;
#else
integer i_indx(char *a, char *b, ftnlen la, ftnlen lb)
#endif
{
    ftnlen i, n;
    char *s, *t, *bend;

    n = la - lb + 1;
    bend = b + lb;

    for(i = 0; i < n; ++i)
    {
        s = a + i;
        t = b;
        while(t < bend)
            if(*s++ != *t++)
                goto no;
        return (i + 1);
    no:;
    }
    return (0);
}
#ifdef __cplusplus
}
#endif
