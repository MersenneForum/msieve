#ifndef DIMODS_H 
#define DIMODS_H


#include <math.h>	/* floor function */
#include "nfstypes.h"

static long      dimods  ARGS((DBLINTC, longc));
                                       /* Nonnegative remainder DBLINT/single */


static long dimods(num, den) DBLINTC num; longc den;
/*
        Divide the (possibly negative) num by positive den,
        returning nonnegative remainder.  Quotient may exceed 2^32.
*/
{
#if DBLINT_REPRESENTATION == DBLINT_LLONG
  register longc rem = num % den;
  return (rem < 0) ? rem + den : rem;
#endif

#if DBLINT_REPRESENTATION == DBLINT_FLOATING
  register DBLINTC quot = floor(num/(DBLINT)den);
/*  printf("dimods: %f %d %f %d\n", num, den, quot); */
  return (long)(num - quot*(DBLINT)den);
#endif
}



#endif /* DIMODS_H */
