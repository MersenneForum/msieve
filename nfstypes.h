#ifndef NFSTYPES_H
#define NFSTYPES_H


#ifdef lint
#ifndef __STDC__
#define __STDC__ 1
#endif
#endif

#include <stdio.h>
//#include <stdlib.h>
//#include <time.h>
#include <math.h>
#include "archtype.h"

/*
	nfstypes.h.  Common definitions for NFS-related programs,
	but excluding declarations dependent on choice of
	multiple precision arithmetic package.
	Suports both  K & R  C  and   ANSI C.
*/

//#include "rcsid.h"
//GENERATE_RCSID(nfstypes_rcsid, "$Id: nfstypes.h,v 1.6 1994/06/07 18:18:56 pmontgom Exp $");

#ifdef __STDC__

#define ARGS(list) list
			/* For function prototypes.			      */
			/* Always insert extra pair of parens.		      */
			/* Example: int myfunc ARGS((intc, double*))	      */
#define NOARGS void
                        /* For functions with no arguments                    */

#else
			/* pre-ANSI C */
#define ARGS(list) ()
#define const
#define signed
#define volatile
			/* Ignore these keywords in pre-ANSI C. */
#define NOARGS

#ifndef labs
#define labs(x) (((x) < 0) ? -(x) : (x))
#endif

#endif  /* __STDC__ */

#ifdef lint
#define do_once(statements) do {statements} while (sqrt(2.0) == 0)
#else 
#define do_once(statements) do {statements} while (FALSE)
#endif

#ifndef MAX
#define MAX(n1, n2) ((n1) < (n2) ? (n2) : (n1))
#define MIN(n1, n2) ((n1) < (n2) ? (n1) : (n2))
					  /* Larger and smaller of n1 and n2. */
	   /* Arguments may be evaluated multiple times - avoid side effects. */
#endif /* MAX */
#ifndef MAX3
#define MAX3(n1, n2, n3) MAX(n1, MAX(n2, n3))
#define MIN3(n1, n2, n3) MIN(n1, MIN(n2, n3))
#endif /* MAX3 */

#define LABS(n) ((n) < 0 ? -(n) : (n))

typedef unsigned char uchar;
/*
	Some versions of the operating system define
	boolean_t, ulong, ushort
	et al in /usr/include/sys/types.h.
	We #include that file now if not previously
	included.  Then we change the names for
	subsequent definitions.
*/

#if ARCH_SGI
#include <sys/types.h>
#define boolean_t MY_BOOLEAN_T
#define     ulong MY_ULONG
#define    ushort MY_USHORT
#endif

typedef unsigned long ulong;
typedef unsigned short ushort;

typedef const uchar  ucharc;
typedef const char   charc;
typedef const int    intc;
typedef const short  shortc;
typedef const long   longc;
typedef const ushort ushortc;
typedef const ulong  ulongc;
typedef const double doublec;
#define CONSTP(type) type*const
		       /* Constant pointer to (possibly non-const) data       */

/* A trailing "c" on a typedef name denotes const data (ANSI C only).	      */
/* For example, charc* ptr might point to a read-only string;		      */
/* the code may modify the subsequently pointer but not the data.	      */
/* Whereas CONSTP(char) ptr might always point to the start of an	      */
/* a string: program can modify data but not pointertring.		      */

typedef signed long psub_t;	/* Index into table of primes		      */
typedef signed long prime_t;	/* Prime for sieving, or value modulo a prime */
typedef psub_t      rowcol_t;	/* Row or column index */

typedef const psub_t psub_tc;
typedef const prime_t prime_tc;
typedef const rowcol_t rowcol_tc;

#if ARCH_CRAY
typedef       int  boolean_t;
				/* No byte-addressing on Cray */
#else
typedef       char boolean_t;	/* TRUE or FALSE value - should be 0 or 1     */
#endif

typedef const char boolean_tc;

#define BOOLEAN boolean_t
#define BOOLEANC boolean_tc
				/* Alternate names for the types */
#define  TRUE ((BOOLEAN)1)
#define FALSE ((BOOLEAN)0)
	/* Problem: Alint complains "Constant in conditional context"
	   for loops beginning "while (TRUE)" and for "do ... while(FALSE)".
	*/
#define BOOL_EQV(b1, b2) ((b1)^(b2)^1)
	/* Return TRUE if b1 and b2 are both TRUE or both FALSE, else FALSE */


typedef struct {
        prime_t p;              /* Norm of ideal (a prime) */
        prime_t id;             /* a/b quotient mod p,
                                   possibly with polynomial
                                   index in upper 4 bits */
                                /* For projective ideal, quotient is p */
        } ideal_t;

#define PID_POLY(pid) ((int)((pid).id >> 28) & 0xf)
                                /* Extract polynomial index from ideal_t pid */
#define PID_IDQUOT(pid) ((pid).id & 0x0fffffff)
				/* Extract a/b mod p from ideal_t pid */
typedef const ideal_t ideal_tc;

typedef struct {
        ideal_t pid;
        int     exponent;
        } ideal_exp_t;

typedef const ideal_exp_t ideal_exp_tc;


/*
	Type DBLINT is intended for signed integers which
	may exceed 2^32 but not 2^53 in absolute value.
	Depending upon implementation,
	it may be a floating type (double) or an integer type (long long).
*/
#define DBLINT_FLOATING 1
#define DBLINT_LLONG 2

#ifndef DBLINT_REPRESENTATION
#define DBLINT_REPRESENTATION DBLINT_FLOATING
#endif

#if DBLINT_REPRESENTATION == DBLINT_FLOATING
typedef double DBLINT;
#endif

#if DBLINT_REPRESENTATION == DBLINT_LLONG
typedef signed long long DBLINT;
#endif

typedef const DBLINT DBLINTC;


typedef char DISTR[20];		/* String large enough for DBLINT */
#define DIPROD(x, y) ((DBLINT)(x)*(DBLINT)(y))

/* See digcd, dimods, diout, each on its own file. */

#if DBLINT_REPRESENTATION == DBLINT_LLONG
#define FITS_IN_DBLINT (fabs(d) < 9223372036854775808.0)	/* 2^63 */
#endif

#if DBLINT_REPRESENTATION == DBLINT_FLOATING
#define FITS_IN_DBLINT(d) (fabs(d) < 9007199254740992.0)	/* 2^53 */
#endif


static int       LOG2    ARGS((ulongc));        /* Truncated base 2 logarithm */
static long      sgcd    ARGS((longc, longc));		       /* Integer GCD */
static boolean_t relprim ARGS((long, long));     /* Test for relatively prime */



static int LOG2(u)
                        /* Base 2 logarithm.  LOG2(0) = -1 */
ulongc u;
{
  static const signed char lowlogs[16] = {-1, 0, 1, 1, 2, 2, 2, 2,
                                           3, 3, 3, 3, 3, 3, 3, 3};
                                          /* Logs of 0 to 15 */
  register ulong u1 = u;
  register int lg2 = 0;
  while (u1 >= 65536) {u1 >>= 16; lg2 += 16;};
  if (u1 >= 256) {u1 >>= 8; lg2 += 8;}
  if (u1 >=  16) {u1 >>= 4; lg2 += 4;}
  return lg2 + (int)lowlogs[u1];
} /* LOG2 */



static long sgcd (m,n)
  longc  m, n;
{
  register long u = labs(m), v = labs(n);
  int ipw2;
  long uorv = u | v;
  if (u == 0) return v;
  ipw2 = 0;
  if ( (uorv & 1) == 0) {
    do {
      ipw2++;
      uorv >>= 1;
    } while ((uorv & 1) == 0);

    u >>= ipw2;
    v >>= ipw2;
  }

  while ((u & 1) == 0) u >>= 1;  /* make u odd */

  while (v != 0){	/* u is odd throughout this loop */
    while ((v&1) == 0) v >>= 1;  /* make v odd */

    while (u > v) {
      u -= v;
      while ((u & 1) == 0) u >>= 1;
    }
    v -= u;
  }

  return (u << ipw2);
}

static boolean_t relprim(u,v)
 register long u, v;
{
  register long t;
  if (!(u&1)) {
    if (!(v&1))  return (0);
    else{  t = u; u = v; v = t; }}

  while (v){	/* u is odd throughout this loop */
    while ((v&1) == 0) v >>= 1;  /* make v odd */
    if ((t=u-v)>0){
      u = v;  v = t;}
    else v -= u;}
  return( u == 1 );
}

#endif /* NFSTYPES_H */
