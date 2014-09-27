
/* Here is basic.c.  Someone who knows gcc can probably make better versions	*/
/* of zdiv21 and zsmulmod (like zaddmulp contribution); be careful 		*/
/* because zsmulmod is also referenced in the new mpqs.c.  			*/

/*
   Numbers are represented in an array of longs.  The format is

   array[0] = array[ZLNGWRD]    length, excluding this word (called lng below)

   array[1] + array[2]*RADIX
            + array[3]*RADIX^2 + ...
            + array[lng]*RADIX^(lng-1)	Absolute value of number.


   Both positive and negaitve numbers are supported.
   A zero has length 0.  Negative numbers have negative lengths.
   Nonzero numbers are normalized, with array[abs(lng)] != 0.

*/


#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "basic.h"

#include "basic2.c"		/* Additional subroutines */

#define LENGTH_ZERO 0

#define	PRIMBND		5000 /* Size of tables used for prime generation */

static basic_t sodinv   ARGS((basic_tc, basic_tc));

static long zaddpos ARGS((MPINPUT, longc, MPINPUT, longc, MPOUTPUT));
			/* Add absolute values */
			/* Assume first length >= second length */

static long zsubpos ARGS((MPINPUT, longc, MPINPUT, longc, MPOUTPUT));
			/* Subtract absolute values, nonnegative difference */

static void zsubmul ARGS((basic_tc, MPMODIFIED, MPINPUT));




/* a = a+t+b*d % RADIX; t = b*d / RADIX */
/* should be optimized */

#if 0
#define		zaddmulp(_a,_b,_d,_t) \
{ register long lb = (_b), ld = (_d), \
		b1 = (_b) & RADIXROOTM, d1 = (_d) & RADIXROOTM; \
  register long aa = *(_a) + b1 * d1; \
  b1 = b1* ( ld >>= NBITSH ) + d1* ( lb >>= NBITSH ) + (aa >> NBITSH); \
  aa = (aa & RADIXROOTM) + ((b1 & RADIXROOTM) << NBITSH) + *(_t); \
  *(_t) = ld*lb + (b1 >> NBITSH) + (aa >> NBITS); \
  *(_a) = (aa & RADIXM); \
}
#endif
#if 1
/*
	This definition of zaddmulp presumes a two's complement
	machine in which integer overflow is ignored
	and where double precision arithmetic is fast.
	The 0.25 allows round to nearest or round towards zero
	(value being rounded should be integer except for roundoff.)
*/

static doublec fradix_inv = 1.0/RADIX;	/* Global constant */



#define zaddmulp(_a,_b,_d,_t) \
{ register basic_tc at = *(_a) + *(_t); \
  register basic_tc aa = (at + (_b)*(_d)) & RADIXM; \
  *(_t) = (long)(0.25 + fradix_inv*(  ((double)(at-aa)) \
		                    + ((double)(_b)) * ((double)(_d)) ) ); \
  *(_a) = aa; \
}
#endif

/* a += d*b, optimize this one */
/* a and b not at overlapping addresses */
/* (except for d=1) */

#define	zaddmul(ams,ama,amb, amblng) \
{ register long lami; \
  register basic_t lams = (ams);  /* Multiplier */ \
  if (lams == 1) { /* Special case in zdiv_help */ \
     lams = 0; \
     for (lami = 1; lami <= amblng; lami++) \
	{  lams += (ama)[lami - 1] + (amb)[lami]; \
	   (ama)[lami - 1] = lams & RADIXM; \
	   lams >>= NBITS;  \
        } \
    (ama)[amblng] += lams; \
  } else { long lamcarry = 0; \
     for (lami = 1; lami <= amblng; lami++) \
	{ zaddmulp(&(ama)[lami - 1], (amb)[lami], lams, &lamcarry);} \
     (ama)[amblng] += lamcarry; \
  } \
}
/* Be careful, the last lama is unnormalized */


/* That was somewhat dirty. If you prefer you can replace the */
/* above macro's (zaddmulp, zaddmul) by the */
/* following C-routines */

/*
zaddmulp(a,b,d,t)
long *a,b,d,*t;
{ register long b1 = b & RADIXROOTM, d1 = d & RADIXROOTM;
  register long aa = *a + b1 * d1;
  b1 = b1* ( d >>= NBITSH ) + d1* ( b >>= NBITSH ) + (aa >> NBITSH);
  aa = (aa & RADIXROOTM) + ((b1 & RADIXROOTM) << NBITSH) + *t;
  *t = d*b + (b1 >> NBITSH) + (aa >> NBITS);
  *a = (aa & RADIXM);
}


zaddmul(d,a,b, blng)
long d,*a,*b, blng;
{ register long i;
  long carry = 0;
  for (i=1; i <= blng; i++)) zaddmulp(&a[i-1], b[i], d, &carry);
  a[blng-1] += carry;
}

*/


/* global variables */

#define IEEE_LENGTH 53
static double IEEE_MAX;		/* 2.0^53 */
static double IEEE_MAXsq;	/* 2.0^106*/
static doublec fradix = (double)RADIX;
static doublec one_plus = 1.0 + 0.5/(double)RADIX;

static double epsilon, fudge;			/* for long division */
static long zseed[ZSIZEP], zranp[ZSIZEP], zprroot[ZSIZEP];	/* for random */

static boolean_t xxprimes[PRIMBND], interval[PRIMBND];	/* to generate primes */
		/* xxprimes[i] is TRUE iff 2*i + 3 is prime */
		/* interval[i] is TRUE iff 2*i + pshift + 3 is prime */
		/* TBD - Merge arrays using multiple bits in one word */

static long pshift = -2;		/* -2 means tables uninitialized */
static long pindex;			/* Location of previous prime    */

/* end global variables */

void zstart()		/* set values global variables */
{ doublec one = (double)1, half = 1/(double)2;
  volatile double fudge1;

  epsilon = one;
#if IEEE_LENGTH == 53
  IEEE_MAX = 1048576.0 * 1048576.0 * 8192.0;	/* 2.0 ^ 53 */
#else
  #error "Compute 2.0^IEEE_LENGTH"
#endif
  IEEE_MAXsq = IEEE_MAX * IEEE_MAX;
  fudge1 = one + epsilon; /* this is sick, but we have to */
  while ( one != fudge1 ) { /* because the test one != (one+epsilon) */
	epsilon = epsilon * half; /* leads to disasters on some machines */
	fudge1 = one + epsilon;
	}
  epsilon += epsilon;
  if ((epsilon*RADIX) > one) printf("decrease RADIX and recompile\n");
  epsilon *= 3;
  fudge = fradix + epsilon*fradix;
  zseed[0] = -1;
}

static void zsubmul(r,a,b)     /* a[s:] -= d*b, optimize this one */
			       /* a and b not at overlapping addresses */
basic_tc r, *b;  basic_t *a;
{ register long rd = RADIX-r, i;
  long carry = RADIX;

  assert (b[ZLNGWRD] > 0);

for (i=(*b++);i>0;i--)
	{ zaddmulp(a,*b,rd,&carry); a++;
	  carry += RADIXM-(*b++);
	}
*a += carry - RADIX; /* unnormalized */
}


void zhalt(e)	/* error handling routine */
		/* See #include file basic.h for error codes */
longc e;
{
  fflush(stdout);
  fprintf(stderr, "zhalt -- fatal error number %ld occurred\n", e);
  fflush(stderr);
  /* Call your favorite error routine here */
  exit(1);
/*NOTREACHED*/  /* to shut up lint */
}


void zcopy(a, b)         /* b = a.  a and b may be the same */
basic_tc *a; basic_t *b;
{ register long i;
  longc lng = LABS(a[ZLNGWRD]);
  for (i = lng; i >= 0; i--) b[i] = a[i];
}


void zcopyabs(a, b)	/* b = |a| .  a and b may be the same*/
basic_tc *a; basic_t *b;
{ register long i;
  longc lng = LABS(a[ZLNGWRD]);
  b[ZLNGWRD] = lng;
  for (i = lng; i > 0; i--) b[i] = a[i];
}





void zintoz(d,a)	/* a = d, d long */
basic_tc d; basic_t a[];
{
  if (d == 0) {
    SET_ZERO(a);
  } else if (d > 0) {
    register ubasic_t dd;
    register long lng;

    lng = 0;
    dd = d;
    do {
      lng++;
      a[lng] = (basic_t)(dd & RADIXM);
      dd >>= NBITS;
    } while (dd != 0);
    a[ZLNGWRD] = lng;
  } else {	/* d < 0 */
    register ubasic_t dd;
    register long lng;

    lng = 0;
    dd = -(ubasic_t)d;
    do {
      lng++;
      a[lng] = (basic_t)(dd & RADIXM);
      dd >>= NBITS;
    } while (dd != 0);
    a[ZLNGWRD] = -lng;
  }
} /* zintoz */


long ztoint(a)	/* returns (long)a, no overflow check */
basic_tc *a;
{ register long d;
  register longc lng = a[ZLNGWRD];

  if (lng == 0) {
    d = 0;
  } else if (lng > 0) {
    register long i;
    d = a[lng];
    for (i=lng-1; i>0; i--) { d = (d << NBITS) + a[i]; }
  } else {	/* lng < 0 */
    register long i;
    d = a[-lng];
    for (i=-lng-1; i>0; i--) { d = (d << NBITS) + a[i]; }
    d = -d;
  }
  return d;
} /* ztoint */


basic_t zcompabs(a, b)	/* |a|>|b|:1; |a|=|b|:0; |a|<|b|:-1 */
basic_tc *a, *b;
{ register longc lnga = LABS(a[ZLNGWRD]);
  register long diff;
  register long i;

  diff = lnga - LABS(b[ZLNGWRD]);
  for (i = lnga; i > 0 && diff == 0; i--) {
    diff = a[i] - b[i];
  }
  return (diff > 0 ? 1 : (diff < 0 ? -1 : 0));
} /* zcompabs */


basic_t zcompare(a,b)	/* a>b:1; a=b:0; a<b:-1 */
basic_tc *a, *b;
{ register longc lnga = a[ZLNGWRD];
  register long diff;
  register long i;

  diff = lnga - b[ZLNGWRD];
  if (lnga >= 0) { /* a >= 0 */
    for (i = lnga; i > 0 && diff == 0; i--) {
      diff = a[i] - b[i];
    }
  } else {	/*  a < 0 */
    register long i;
    for (i = -lnga; i > 0 && diff == 0; i--) {
      diff = b[i] - a[i];
    }
  } /* if lnga */
  return (diff > 0 ? 1 : (diff < 0 ? -1 : 0));
} /* zcompare */


void zsadd(a, d, b) /* b = a+d, abs(d)<RADIX, output can be input */
basic_tc a[],d; basic_t b[];
{
  basic_t x[2];
  x[ZLNGWRD] = 1;
  if (d > 0) {
    x[1] = d;
    zadd(a, x, b);
  } else if (d < 0) {
    x[1] = -d;
    zsub(a, x, b);
  } else {
    zcopy(a, b);
  }
}



static long zaddpos(a, lnga, b, lngb, c)
longc lnga, lngb;
basic_tc *a, *b;
basic_t *c;
{	/* Add absolute values, return length of sum */
	/* Assume longer operand is in a */
  register long i;
  register ubasic_t carry;
  if (lnga < lngb || lngb < 0) zhalt(ZHALT_ZADDPOS_BAD_LENGTHS);

  carry = 0;
  for (i = 1; i <= lngb; i++) {
    carry += a[i] + b[i];
    c[i] = carry & RADIXM;
    carry >>= NBITS;
  }

  for (i = lngb + 1; i <= lnga; i++) {
    carry += a[i];
    c[i] = carry & RADIXM;
    carry >>= NBITS;
  }

  if (carry  == 0) {
    return lnga;
  } else {
    c[lnga + 1] = carry;
    return lnga+1;
  }
} /* zaddpos */



static long zsubpos(a, lnga, b, lngb, c)
longc lnga, lngb;
basic_tc *a, *b;
basic_t *c;
{	/* Subtract absolute values, return length of difference */
	/* Difference assumed to be nonnegative. */
  register long i;
  register ubasic_t carry;
  if (lnga < lngb || lngb < 0) zhalt(ZHALT_ZSUBPOS_BAD_LENGTHS);

  carry = 1;			/* 0 if borrowing, else +1 */
  for (i = 1; i <= lngb; i++) {
    carry += (ubasic_t)(RADIXM + a[i] - b[i]);
    c[i] = carry & RADIXM;
    carry >>= NBITS;
  }

  for (i = lngb + 1; i <= lnga; i++) {
    carry += (RADIXM + a[i]);
    c[i] = carry & RADIXM;
    carry >>= NBITS;
  }

  if (carry == 0) {
    zhalt(ZHALT_ZSUBPOS_NEGATIVE);
  }

  i = lnga;
  while (i > LENGTH_ZERO && c[i] == 0) i--;
  return i;
} /* zsubpos */



void zadd(a, b, c)	/* c = a+b; output can be input */
basic_tc *a,*b; basic_t c[];
{ register long lnga, lngb;
  register basic_tc *pa, *pb;

  if (a[ZLNGWRD] >= b[ZLNGWRD]) {	/* Larger LNGWRD -> pa */
    pa = a; lnga = a[ZLNGWRD];
    pb = b; lngb = b[ZLNGWRD];
  } else {
    pa = b; lnga = b[ZLNGWRD];
    pb = a; lngb = a[ZLNGWRD];
  }

  if (lngb >= 0) {		/* Both non-negative */
    c[ZLNGWRD] = zaddpos(pa, lnga, pb, lngb, c);
  } else if (lnga <= 0) {	/* Both non-positive */
    c[ZLNGWRD] = -zaddpos(pb, -lngb, pa, -lnga, c);  /* -(-b + -a) */
  } else if (zcompabs(pa, pb) >= 0) {  /* a >= -b > 0 */
    c[ZLNGWRD] = zsubpos(pa, lnga, pb, -lngb, c);
  } else {			       /* -b > a > 0 */
    c[ZLNGWRD] = -zsubpos(pb, -lngb, pa, lnga, c);
  }
}

void zsub(a, b, c)	/* c = a-b; output can be input */
			/* Same as zadd except negated b[ZLNGWRD] */
basic_tc *a,*b; basic_t c[];
{ register long lnga, lngb;
  register basic_tc *pa, *pb;

  if (a[ZLNGWRD] + b[ZLNGWRD] >= 0) {
    pa = a; lnga = a[ZLNGWRD];
    pb = b; lngb = -b[ZLNGWRD];
  } else {
    pa = b; lnga = -b[ZLNGWRD];
    pb = a; lngb = a[ZLNGWRD];
  }

  if (lngb >= 0) {
    c[ZLNGWRD] = zaddpos(pa, lnga, pb, lngb, c);
  } else if (lnga <= 0) {
    c[ZLNGWRD] = -zaddpos(pb, -lngb, pa, -lnga, c);
  } else if (zcompabs(pa, pb) >= 0) {
    c[ZLNGWRD] = zsubpos(pa, lnga, pb, -lngb, c);
  } else {
    c[ZLNGWRD] = -zsubpos(pb, -lngb, pa, lnga, c);
  }
}


void zsmul(a, d, b)  /* b = d*a, d < RADIX; output can be input */
basic_tc *a,d; basic_t b[];
{ 
  long slnga = a[ZLNGWRD];
  if (d == 0 || IS_ZERO(a) ) {
    SET_ZERO(b);
  } else { 
    register long i, lngb;
    long  dabs;
    dabs = d;
    if (d < 0) {
      dabs = -d;
      slnga = -slnga;
    }
    lngb = LABS(slnga); 

    for (i = 1; i <= lngb; i++) b[i] = a[i];
			/* Caution -- a and b may be same array */
    b[lngb+1] = 0;
    
    zaddmul(dabs-1, &b[1], &b[0], lngb);

    lngb++; 
    while (b[lngb] == 0) lngb--; 
    assert (lngb > 0);

    b[ZLNGWRD] = (slnga >= 0 ? lngb : -lngb);
  } /* d != 0 */
} /* zsmul */


void zmul(a, b, c)	/* c = a*b; output cannot be input */
basic_tc *a,*b; basic_t c[];
{ register long lnga = a[ZLNGWRD], lngb = b[ZLNGWRD];
  if (IS_ZERO(a) || IS_ZERO(b)) {
      SET_ZERO(c);
  } else {
      boolean_t neg_product;
      register basic_tc *pa, *pb;
      register long i, lngc;

      if (lngb < 0) { lnga = -lnga; lngb = -lngb;}
      neg_product = lnga < 0;
      if (neg_product) lnga = -lnga;

      lngc = lnga + lngb;
      for (i = 1; i <= lngc; i++) c[i] = 0;

      if (lnga >= lngb) {	/* Put longer operand in pa */
	pa = a; pb = b;
      } else {
	pa = b; pb = a; lnga = lngc - lnga; lngb = lngc - lngb;
      }
      for (i = 1; i <= lngb ; i++) {
	zaddmul(pb[i], &c[i], pa, lnga);
      }

      while (c[lngc] == 0) lngc--;
      c[ZLNGWRD] = (neg_product ? -lngc : lngc);
  }
}

/*
	zdiv21 returns quot, rem so

	quot = (numhigh*RADIX + numlow)/denom;
	rem  = (numhigh*RADIX + numlow)%denom;

Assumes 0 <= numhigh < denom < RADIX and 0 <= numlow < RADIX.
*/

#define zdiv21(numhigh, numlow, denom, quot, rem) {\
    register long lr21, lq21 = (long)(((fradix*(double)(numhigh)) \
		              + (double)(numlow)) / (double)(denom)); \
    if (1) { \
      long lprodhigh = 0, lprodlow = 0; \
      zaddmulp(&lprodlow, lq21, denom, &lprodhigh); \
      lr21 = RADIX*(numhigh - lprodhigh) + (numlow - lprodlow); \
    } else { \
/* Following works in many two's complement architectures. */ \
      lr21 = (numhigh << NBITS) + numlow  - lq21*denom; \
    }\
    if (lr21 < 0) { do {lq21--;} while ((lr21 += denom) < 0); \
    } else      { while (lr21 >= denom) {lr21 -= denom; lq21++;}; \
    }  \
    *quot = lq21; *rem = lr21; \
}



basic_t zsdiv(a, d, b)	/* b = a/d, returns a%ld; output can be input */
basic_tc a[], d; basic_t b[];
{
    longc slnga = a[ZLNGWRD]; 
    register long lnga = LABS(slnga);
    register long den = d;
    if (IS_ZERO(a)) {	/* Numerator is zero */
        SET_ZERO(b);
        return 0;
    } else if (den >= RADIX || den < 0) {
        long zd[ZSIZEP];
        zintoz(den, zd); 
	zdiv(a, zd, b, zd); 
        return ztoint(zd);
    } else if (lnga == 1) {
        longc numer = (slnga >= 0 ? a[1] : -a[1]);
        basic_t q = numer/den;
        basic_t carry = numer - q*den;
        if (carry < 0) {
	  q--;
          carry += den;
	}
        zintoz(q, b);
        return carry;
    } else {
        register long i, carry;
        carry = a[lnga];
        if (carry < den) {
            lnga--;
        } else {
            carry = 0;
        }
        if (den < (1L << (IEEE_LENGTH - NBITS))) {	
				/* d*RADIX fits in a double */
            register doublec fd = (double)den;
	    register doublec fdinv = one_plus/fd;
	    register double  fcarry = (double)carry;

	    for (i = lnga; i > 0; i--) {
		register basic_t iquot;
	        fcarry = fcarry*fradix + (double)a[i];
		iquot = (basic_t)(fcarry*fdinv);
				/* Exact or one too big */
	        fcarry -= fd*(double)iquot;
	        b[i] = iquot;
                if (fcarry < 0.0) {
		    b[i]--;
		    fcarry += fd;
		}
	    }
	    carry = (basic_t)fcarry;
	} else if (den < RADIXROOT) {
	    register longc rdivd = RADIX/den, rmodd = RADIX - den*rdivd;
	    for (i = lnga; i > 0; i--)
		{ register longc temp = a[i] + rmodd*carry;
		b[i] = rdivd*carry + temp/den;
		carry = temp % den;
		}
	} else {
	    for (i = lnga; i > 0; i--)
	        { long newcarry;
		zdiv21(carry, a[i], den, &b[i], &newcarry);
		carry = newcarry;
		}
        }

	while (b[lnga] == 0) lnga--;	/* (Quotient cannot be zero) */

        assert(lnga > 0);

        if (slnga < 0) {		/* Numerator was negative */
           b[ZLNGWRD] = -lnga;
           if (carry != 0) {	/* Ensure remainder nonnegative */
             carry = den - carry;
	     zsadd(b, -1L, b);
	   }
	} else {
          b[ZLNGWRD] = lnga; 
	}
	return carry ;
    } /* else */
} /* zsdiv */



void zread(a)		/* read a (nonnegative) from standard input */
long a[];
{ char in[ZDECSIZE];
  register long d = 0;
  scanf("%s",in);
  a[0] = 1; a[1] = 0;
  while (in[d] != '\0')
	{ if (in[d] == '\\')
		{ scanf("%s",in); d = 0; }
	  else
		{ zsmul(a,(long)10,a);
		  zsadd(a,(long)(in[d++]-'0'),a);
		};
}	}

void fzread(fil,a)		/* read a (nonnegative) from file fil */
FILE *fil;
long a[];
{ char in[ZDECSIZE];
  register long d = 0;
  fscanf(fil,"%s",in);
  a[0] = 1; a[1] = 0;
  while (in[d] != '\0')
	{ if (in[d] == '\\')
		{ fscanf(fil,"%s",in); d = 0; }
	  else
		{ zsmul(a,(long)10,a);
		  zsadd(a,(long)(in[d++]-'0'),a);
		};
}	}

long zwrite(a)	/* write a on standard output */
basic_tc a[];
{ long ca[ZZSIZEP], out[ZDECSIZE];
  register long sa, i;
  long result;
  if (a[ZLNGWRD] < 0) printf("-");
  zcopyabs(a, ca);
  i = -1;
  do out[++i] = zsdiv(ca,(long)10,ca) ;
  while (ca[0] > 1 || ca[1] != 0 );
  sa = 0;
  if (i < 0 ) i = 0;
  result = i+1;
  while (i >= 0 )
	{ printf("%ld",out[i--]);
	  sa ++;
	  if ((i >= 0) && (sa == 70 ))
		{ sa = 0; printf("\\\n"); }
	}
  return(result);
}

long fzwrite(fil,a,lll,eol,bol)
		/* write a on file fil, at most lll digits per line */
		/* end line with eol, begin next line with bol */
FILE *fil;
basic_tc a[],lll;
charc *eol, *bol;
{ long ca[ZZSIZEP], out[ZDECSIZE];
  register long sa, i;
  long result;
  if (a[ZLNGWRD] < 0) fprintf(fil, "-");
  zcopyabs(a, ca);
  i = -1;
  do { 
    out[++i] = zsdiv(ca,(long)10,ca) ;
  } while (ca[0] > 1 || ca[1] != 0 );
  sa = 0;
  if (i < 0 ) i = 0;
  result = i+1;
  while (i >= 0 )
	{ fprintf(fil,"%ld",out[i--]);
	  sa ++;
	  if ((i >= 0) && (sa == lll ))
		{ sa = 0; fprintf(fil,"%s\n%s",eol,bol); }
	}
  return(result);
}



void zdiv_help(a, lnga, b, lngb, q, r)
	/* q = a/b, r = a%b, a and b positive */
	/* Arrays c, d should be one longer than larger of a, b */
basic_tc a[], b[]; basic_t q[], r[];
longc lnga, lngb;
{
  basic_tc btop = b[lngb];
  register long i;

  for (i = 1; i <= lnga; i++) r[i] = a[i];
  r[ZLNGWRD] = lnga;

  if (lngb == 1 && btop < RADIX) {
     zintoz(zsdiv(r, btop, q), r);
  } else if (lnga < lngb) {		/* Implies quotient is zero */
     SET_ZERO(q);
  } else {
	  long sq = lnga - lngb + 1, sr;
          double btopinv;
          /* Dividend in r -- indices 0 to lnga + 1 = sq + lngb  */
	  /* Quotient in q -- indices 0 to sq                */

          btopinv = fradix*(double)btop;
	  if (lngb > 1) btopinv += (double)b[lngb-1];
	  btopinv *= fradix;
	  if (lngb > 2) btopinv += (double)b[lngb-2];
	  btopinv = fudge / btopinv;

          r[lnga + 1] = 0;

	  for (i = sq; i > 0; i--) {
	      CONSTP(basic_t) pr = &r[i + lngb];

	      if (pr[0] >= 0 || pr[-1] >= btop) {
	          register long qq;
	          double aux;

		  aux = fradix*(double)pr[0] + (double)pr[-1];
                  aux = 1.0 + fradix*aux;
		  if (pr >= &r[2]) aux += (double)pr[-2];
		  qq = (long)(btopinv*aux + 0.5);

		/* dirty, but safe. On most machines (long)(btopinv*aux) */
		/* is correct, or one too big; on some however it becomes*/
		/* too small. Could change zstart, but +0.5 and a while */
		/* instead of one single IF is safer */

	          zsubmul(qq,&r[i],&b[0]);
		  while (pr[0] < 0) {	/* possibly unnormalized */
		      qq--;
		      zaddmul(1, &r[i], &b[0], b[0]);
	          }
		  q[i] = qq;
	      } else {
		  q[i] = 0;
	      }
	  } /* loop on i */

          sr = lngb;
	  while (sr > LENGTH_ZERO && r[sr] == 0) sr--;
	  r[0] = sr;

	  while (sq > LENGTH_ZERO && q[sq] == 0) sq--;
	  q[0] = sq;

	} /* non-trivial case */
}


void zdiv(a, b, q, r)	/* q = a/b, r = a%b, output can be input */
			/* Allow double-length operends */
basic_tc a[], b[]; basic_t q[], r[];
{
  basic_t qtemp[ZZSIZEP + 1], rtemp[ZZSIZEP + 1];
  longc slnga = a[ZLNGWRD], slngb = b[ZLNGWRD];

  if (IS_ZERO(b)) zhalt(ZHALT_ZDIV_ZERO_DIVISOR);
  zdiv_help(a, LABS(slnga), b, LABS(slngb), qtemp, rtemp);

  if (slnga < 0) {
    qtemp[ZLNGWRD] = -qtemp[ZLNGWRD]; 
    rtemp[ZLNGWRD] = -rtemp[ZLNGWRD];
  }
  if (slngb < 0) { 
    qtemp[ZLNGWRD] = -qtemp[ZLNGWRD]; 
  }
  if (IS_NEG(qtemp) && !IS_ZERO(rtemp)) { /* Give remainder same sign as b */
    zadd(rtemp, b, r);
    zsadd(qtemp, -1L, q);
  } else {
    zcopy(rtemp, r);
    zcopy(qtemp, q);
  }
}


void zmod(a, b, r)	/* r = a%b, output can be input */
			/* Remainder has same sign as divisor */
			/* Allow double-length operands */
basic_tc a[], b[]; basic_t r[];
{
  basic_t qtemp[ZZSIZEP + 1], rtemp[ZZSIZEP + 1];
  longc slnga = a[ZLNGWRD], slngb = b[ZLNGWRD];

  if (IS_ZERO(b)) zhalt(ZHALT_ZMOD_ZERO_DIVISOR);
  zdiv_help(a, LABS(slnga), b, LABS(slngb), qtemp, rtemp);

  if (slnga < 0) {
    qtemp[ZLNGWRD] = -qtemp[ZLNGWRD]; 
    rtemp[ZLNGWRD] = -rtemp[ZLNGWRD];
  }
  if (slngb < 0) { 
    qtemp[ZLNGWRD] = -qtemp[ZLNGWRD]; 
  }
  if (IS_NEG(qtemp) && !IS_ZERO(rtemp)) {
    zadd(rtemp, b, r);
  } else {
    zcopy(rtemp, r);
  }
} /* zmod */


void zaddmod(a,b,n,c)	/* c = (a+b)%n, with 0<=a,b<n */
basic_tc a[],b[],n[]; basic_t c[];
{
  zadd(a, b, c);
  if (zcompare(c, n) >= 0 ) zsub(c, n, c);
}


void zsubmod(a,b,n,c)	/* c = (a-b)%n, with 0<=a,b<n */
basic_tc a[],b[],n[]; basic_t c[];
{
  zsub(a, b, c);
  if (IS_NEG(c)) zadd(c, n, c);
}


long zsmulmod(a, b, n)	/* return (a*b)%n */
			/* Assume |a| < n and |b| < n */
longc a,b,n;
{
  register doublec dab = (double)a * (double)b;
  register doublec dn = (double)n;
  register doublec dquot = dab/dn;
		          /* Divide time can overlap other operations */
  register long lr;
  if (dab * dab < IEEE_MAXsq)  {   /* Product was exact */
				/* Avoid fabs, which is not done inline */
    lr = (long)(dab - dn * (double)(long)dquot);
  } else {
/* Many machines compute the following modulo 2^32, which is OK */
    lr = a*b - n*(long)dquot;
  }
  do {lr -= n;} while (lr >= 0);
  do {lr += n;} while (lr < 0);
  return lr;
}


void zmulmod(a,b,n,c)	/* c = (a*b)%n */
basic_tc a[],b[],n[]; basic_t c[];
{ long mem[ZZSIZEP];
  zmul(a,b,mem);	zmod(mem,n,c);
}


void zmstart(n, mont)	/* set up values for Montgomery's multiplication */
basic_tc n[];
CONSTP(zmont_t) mont;
{ basic_t x[ZSIZEP+1], i;
  ubasic_t zminv, zminvrem;
  basic_tc n1 = n[1];
  basic_tc zmtop = n[0];
  if (zmtop >= ZSIZEP || zmtop <= 0) zhalt(ZHALT_ZMSTART_TOO_LONG);
  if ((n1&1) == 0) zhalt(ZHALT_ZMSTART_EVEN_MODULUS);


  zminv = 1;
  zminvrem = (n1 + 1) >> 1;
  for (i = 1; i < NBITS; i++) {
        /* At this point     zminv*n1 = 2^i * zinvrem - 1
	             and   0 < zminv < 2^i.
           Replace i by i+1.
        */
    if ((zminvrem & 1) != 0) {
      zminv += (1L << i);
      zminvrem += n1;
    }
    zminvrem >>= 1;
  }			/* Now zminv*n1 = RADIX*zinvrem - 1 */

  mont->zminv = zminv;		/* -1/n[1] mod RADIX */
  for (i=1; i<=zmtop; i++) x[i] = 0;
  x[zmtop+1] = 1; x[0] = zmtop+1;
  zmod(x, n, mont->one);        /* one   =  RADIX^zmtop mod n */
  zmulmod(mont->one, mont->one, n, mont->tom);
				/* tom   =  RADIX^(2*zmtop) mod n */
  zcopy(n, mont->n);
}

void ztom(a,mont,b)	/* b is Montgomery representation of a */
basic_tc a[]; basic_t b[];
CONSTP(zmont_tc) mont;
{
  if (zcompare(mont->n, a) <= 0 || IS_NEG(a) ) {
    zmod(a, mont->n, b);
    zmont(b, mont->tom, mont, b);
  } else {
    zmont(a, mont->tom, mont, b);
  }
}

void zmtoz(a,mont,b)	/* Montgomery a is backtransformed to normal b */
basic_tc a[];
basic_t b[];
CONSTP(zmont_tc) mont;
{ long x[2];
  x[0] = 1; x[1] = 1;
  zmont(a, x, mont, b);
}

/*
   Macro to add quadwords.  Output can overlap first input but not second.
*/

#define QUAD_ADD(q1, q2, q3) {q3[0] = q1[0] + q2[0];  \
	q3[1] = q1[1] + q2[1] + ((ubasic_t)q3[0] < (ubasic_t)q2[0]);}

#define QUAD_UPPER(q) ((q[1] << (32 - NBITS)) | ((ubasic_t)q[0] >> NBITS))

#define QUAD_ADD_SCALAR(q, s) {q[0] += (ubasic_t)(s);  \
		q[1] += ((ubasic_t)q[0] < (ubasic_t)(s));}

#define QUAD_ZERO(q) q[0] = q[1] = 0
#define QUAD_SET_SCALAR(q, s) {q[0] = (ubasic_t)(s); q[1] = 0;}


void zmont(a, b, mont, c)  /* c is mont product of mont's a and b */
		           /* output can be input */
basic_tc a[], b[];
CONSTP(zmont_tc) mont;
basic_t c[];
{
  register int i, j;
  CONSTP(basic_tc) pa = (b[0] < a[0] ? a : b);	/* Longer operand */
  CONSTP(basic_tc) pb = (b[0] < a[0] ? b : a);  /* Shorter operand */
  intc alng = pa[0], blng = pb[0], nlng = mont->n[0];
  ubasic_t acopy[ZSIZEP], tmp[ZSIZEP];
  register ubasic_t carry, mula, muln;
  ubasic_t quad[2];

  mula = pb[1];
  QUAD_SET_SCALAR(quad, 0);
  acopy[1] = pa[1];
  PRDPPA(mula, acopy[1], quad);
  muln = (quad[0] * mont->zminv) & RADIXM;   /* Bottom bits of product */
  PRDPPA(muln, mont->n[1], quad);
  assert ((quad[0] & RADIXM) == 0);
  carry = QUAD_UPPER(quad);

  if (nlng == 1) {
     tmp[nlng] = carry;
  } else if (pa == pb)  {			/* Squaring */

     for (i = 2; i <= nlng; i++) {
		    /*  tmp = (mula*a[1] + muln*n)/RADIX  */

	acopy[i] = (i <= alng ? pa[i] : 0);
        QUAD_SET_SCALAR(quad, carry);
        PRDPPA(muln, mont->n[i], quad);
        tmp[i-1] = quad[0] & RADIXM;
        carry = QUAD_UPPER(quad);
     }
     tmp[nlng] = carry;


     for (j = 2; j <= nlng; j++) {
	boolean_tc bigmul = (acopy[j] >= RADIXH);
        register ubasic_t carry, mula, muln;
        ubasic_t quad[2];

	mula = (acopy[j] << 1) & RADIXM;
	QUAD_SET_SCALAR(quad, tmp[1]);
	PRDPPA(mula, acopy[1], quad);
        muln = (quad[0] * mont->zminv) & RADIXM;
	PRDPPA(muln, mont->n[1], quad);
	assert ((quad[0] & RADIXM) == 0);
	carry = QUAD_UPPER(quad);

	if (bigmul) {
	   carry += acopy[1];
	   for (i = 2; i < j; i++) {
	      QUAD_SET_SCALAR(quad, carry);
	      QUAD_ADD_SCALAR(quad, tmp[i]);
	      PRDPPA(mula, acopy[i], quad);
	      PRDPPA(muln, mont->n[i], quad);
	      tmp[i-1] = quad[0] & RADIXM;
	      carry = QUAD_UPPER(quad) + acopy[i];
	   }
	} else {		/* !bigmul */
	   for (i = 2; i < j; i++) {
	      QUAD_SET_SCALAR(quad, carry);
	      QUAD_ADD_SCALAR(quad, tmp[i]);
	      PRDPPA(mula, acopy[i], quad);
	      PRDPPA(muln, mont->n[i], quad);
	      tmp[i-1] = quad[0] & RADIXM;
	      carry = QUAD_UPPER(quad);
	   }
	}
	QUAD_SET_SCALAR(quad, carry);
	PRDPPA(acopy[j], acopy[j], quad);
	for (i = j; i <= nlng; i++) {
	   PRDPPA(muln, mont->n[i], quad);
	   QUAD_ADD_SCALAR(quad, tmp[i]);
           tmp[i-1] = quad[0] & RADIXM;
	   quad[0] = QUAD_UPPER(quad);
	   quad[1] = 0;
	}
        tmp[nlng] = quad[0];

     } /* for j */

  } else {				/* General multiplication */
     for (i = 2; i <= nlng; i++) {
		    /* tmp = (mula*a + muln*n)/RADIX */
	basic_tc ia = (i <= alng ? pa[i] : 0);

        acopy[i] = ia;
        QUAD_SET_SCALAR(quad, carry);
        PRDPPA(mula, ia, quad);
        PRDPPA(muln, mont->n[i], quad);
        tmp[i-1] = quad[0] & RADIXM;
        carry = QUAD_UPPER(quad);
     }
     tmp[nlng] = carry;

     for (j = 2; j <= nlng; j++) {
         mula = (j <= blng ? pb[j] : 0);
         QUAD_SET_SCALAR(quad, tmp[1]);
         PRDPPA(mula, acopy[1], quad);
         muln = (quad[0] * mont->zminv) & RADIXM;   /* Bottom bits of product */
         PRDPPA(muln, mont->n[1], quad);
         assert ((quad[0] & RADIXM) == 0);
         carry = QUAD_UPPER(quad);

         for (i = 2; i <= nlng; i++) {
			    /* tmp = (tmp + mula*a + muln*n)/RADIX */
            QUAD_SET_SCALAR(quad, tmp[i] + carry);
	    PRDPPA(mula, acopy[i], quad);
            PRDPPA(muln, mont->n[i], quad);
            tmp[i-1] = quad[0] & RADIXM;
            carry = QUAD_UPPER(quad);
         }
         tmp[nlng] = carry;		    /* Caution - this may overflow */
      }
  }

  for (i = 1; i <= nlng; i++) {
     c[i] = tmp[i];
  }
  i = nlng;
  while (i > LENGTH_ZERO && tmp[i] == 0) i--;	/* Determine new length */
  c[0] = i;
  if (zcompare(c, mont->n) >= 0) zsub(c, mont->n, c);
}

void zexp(a,e,n,b)	/* b = (a^e)&n, b and a can be identical, but */
			/* both unequal to n and e */
			/* only if RADIX is power of two */
basic_tc a[],e[],n[]; basic_t b[];
{ register long i,j, cntr, ei;
  basic_t sq[ZZSIZEP];
  basic_tc *pae;
  if ((e[0] == 1) && (e[1] == 0 ))
	{ b[0] = 1; b[1] = 1; }
  else
	{ zmod(a,n,sq);
	  b[0] = 1; b[1] = 1; cntr = 0; pae = &e[0];
	  for (i=(*pae); i>0; i--)
		{ for (j=cntr; j>0; j--) zmulmod(sq,sq,n,sq);
		  ei = *(++pae); cntr = NBITS;
		  if (ei > 0)
			{ while (ei > 1)
				{ if ( (ei&1) == 1 ) zmulmod(b,sq,n,b);
				  ei >>= 1;
				  zmulmod(sq,sq,n,sq); cntr--;
				}
			  zmulmod(b,sq,n,b);
			}
		}
	}
}

long zsexp(a,e,n)	/* Single precision version of zexp */
longc a,e,n;
{
    if (e == 0)
	return(1);
    else
	{ register long aa = a%n, ipow2 = 1;
	  long b = aa;
	  while ((ipow2 << 1) <= e) ipow2 <<= 1;
	  while ((ipow2 >>= 1) != 0) {
		/* Answer = (b^(2*ipow2))*(aa^(e%(2*ipow2))) mod n */
		b = zsmulmod(b, b, n);
		if ((e & ipow2) != 0) b = zsmulmod(b, aa, n);
		}
	  return(b);
	}
}

/* returns the inverse of n modulo the prime p	*/
/* first the old version, undefined if not coprime */
/* then a faster binary version */
/*
long sinv(n,p)
long n,p;
{ register long q=n, r, nq, nr, u=1, w, nw, par = 0;
if (n > p) { w = 0; r = p; }
else r = p- (w = p/n) *n;
nr = r;
while (nr != 0) {
	if ( (nr = q - r) < r) nq = 1;
	else if ( ( nr -= r ) < r) nq = 2;
	else nr = q - ( nq = q/r ) *r;
	nw = nq*w+u;	u = w;
	w = nw;		q = r;
	r = nr;		par = 1-par;
	}
if (par == 0) return(u);
else return(p-u);
}
*/

static basic_t sodinv(n,p) longc n, p;
{ /* only for odd p>=3 and n>0 */
	register long n1 = n, n2 = p, preg = p, m1 = 1, m2 = 0;
		/* m1*n == n1 mod p,  m2*n == n2 mod p */
                /* n2 is odd, n1 is odd after initialization */
	while (!(n1&1)) {
		n1 >>= 1;
		if (m1&1) m1 += preg;
		m1 >>= 1;
	}
	if (n1 == 1) return m1;
	while (n1 != n2) {
		if (n1 >n2) {
			n1 -= n2; m1 -= m2;
			if (m1 < 0) m1 += preg;
			do {
				n1 >>= 1;
				if (m1&1) m1 += preg;
				m1 >>= 1;
			} while (!(n1 & 1));
			if (n1 == 1) return m1;
		} else {
			n2 -= n1; m2 -= m1;
			if (m2 < 0) m2 += preg;
			do {
				n2 >>= 1;
				if (m2&1) m2 += preg;
				m2 >>= 1;
			} while (!(n2 & 1));
			if (n2 == 1) return m2;
		}
	}
	zhalt(ZHALT_SODINV_NOT_COPRIME); /* non-trivial gcd n1 */
/*NOTREACHED*/	/* to shut up lint */
}


long sinv (n,p) longc n,p; { /* halts if non-trivial gcd or even p */
	if (p < 3 || (!(p & 1)) || n <= 0) zhalt(ZHALT_SINV_INVALID_ARGUMENTS);
	return( sodinv(n,p) );
}

long zinv(ain,nin,inv)	/* if 0, then (ain*inv)%nin=1 */
			/* if 1, then gcd(ain,nin)=inv > 1 */
basic_tc ain[],nin[]; basic_t inv[];
{ basic_t a[ZZSIZEP], n[ZZSIZEP], q[ZZSIZEP], u[ZZSIZEP],
      w[ZZSIZEP], x[ZZSIZEP], y[ZZSIZEP], z[ZZSIZEP],
      diff, ilo, sa, sn, temp,
      try11, try12, try21, try22, got11, got12, got21, got22,
      fast, parity, gotthem, *p;
  double hi, lo, dt, fhi, flo, num, den;


  fhi = 1.0 + epsilon; flo = 1.0 - epsilon;
  zcopyabs(nin, n);
  zmod(ain, n, a);
  SET_ONE(u);
  SET_ZERO(w);
  while (!IS_ZERO(n)) 
	{
	  gotthem = 0; sa = a[0]; sn = n[0]; diff = sa-sn;
	  if (diff == 0 || diff == 1 )
		{ sa = a[0]; p = &a[sa];
		  num = (double)(*p)*fradix;
		  if (sa > 1) num += (*(--p));
		  num *= fradix;
		  if (sa > 2) num += (*(p-1));
		  sn = n[0]; p = &n[sn];
		  den = (double)(*p)*fradix;
		  if (sn > 1) den += (*(--p));
		  den *= fradix;
		  if (sn > 2) den += (*(p-1));
		  hi = fhi*(num+1.0)/den;
		  lo = flo*num/(den+1.0);
		  if (diff > 0) { hi *= fradix; lo *= fradix; }
		  try11 = 1; try12 = 0; try21 = 0; try22 = 1;
		  parity = 1; fast = 1;
		  while (fast > 0)
			{ parity ^= 1;
/* be careful, some mailers might change the previous line into */
/* open curly bracket parity tilde-equal one semicolon */
/* the correct line is */
/* open curly bracket parity hat-equal one semicolon */
			  if (hi >= fradix )
				{ fast = 0; }
			  else	{ ilo = (long)lo;
				  if (ilo == 0 || ilo < (long)hi)
					{ fast = 0; }
				  else	{ dt = lo; lo = flo/(hi-ilo);
					  if (dt > ilo) hi = fhi/(dt-ilo);
					  else hi = fradix;
					  temp = try11; try11 = try21;
					  if ((RADIX-temp)/ilo < try21)
						fast = 0;
					  else try21 = temp+ilo*try21;
					  temp = try12; try12 = try22;
					  if ((RADIX-temp)/ilo < try22)
						fast = 0;
					  else try22 = temp+ilo*try22;
					  if ((fast > 0) && (parity >0))
						{ gotthem = 1;
						  got11 = try11; got12 = try12;
						  got21 = try21; got22 = try22;
						}
					}
				}
			}
		}
	  if (gotthem > 0)
		{ zsmul(u,got11,x); zsmul(w,got12,y);
		  zsmul(u,got21,z); zsmul(w,got22,w);
		  zadd(x,y,u);     zadd(z,w,w);
		  zsmul(a,got11,x); zsmul(n,got12,y);
		  zsmul(a,got21,z); zsmul(n,got22,n);
		  zsub(x,y,a);     zsub(n,z,n);
		}
	  else	{ zdiv(a,n,q,a); zmul(q,w,x); zadd(u,x,u);
		  if (!IS_ZERO(a))
			{ zdiv(n,a,q,n); zmul(q,u,x); zadd(w,x,w); }
		  else	{ 
			  zcopy(n, a);
			  SET_ZERO(n);  
			  zsub(nin,w,u);
			}
		}
	}
  if (IS_ONE(a)) {
    zcopy(u, inv);		/* GCD = 1, return multiplicative inverse */
    return 0; 
  } else {			/* Non-trivial GCD */
    zcopy(a, inv); 
    return 1;
  }
}

void zrstart(s)	/* initalize globals zseed = s, zranp = 2^89-1, */
		/* zprroot = 3^101 mod zranp */
basic_tc s;
{ register long i, ones, res;
  long e[2];
  zintoz(s,zseed);
  ones = 89/NBITS; res = 89%NBITS;
  for (i=1; i<=ones; i++) zranp[i] = RADIXM;
  if (res != 0) zranp[++ones] = ( RADIXM >> (NBITS - res) );
  zranp[0] = ones;
  zintoz((long)3,zprroot); e[0] = 1; e[1] = 101;
  zexp(zprroot,e,zranp,zprroot);
}

long zrandom(b) /* returns random in [0,b) */
longc b;
{ long zb[ZSIZEP];
  if (zseed[0] < 0) zhalt(ZHALT_ZRANDOM_UNINITIALIZED);
  zmulmod(zseed,zprroot,zranp,zseed);
  zintoz(MAX(b, 1) ,zb);
  zmod(zseed,zb,zb);
  return(ztoint(zb));
}

boolean_t zprime(m,t)		/* 0 if m not prime */
			        /* 1 if m passes t pseudo prime tests */
basic_tc m[];
intc t;
{
 if ((m[0] == 1) && (m[1] <= RADIXROOT)) {
	  basic_tc sm = m[1];
          basic_t i;
	  if (sm < 2) return FALSE;
	  if (sm == 2) return TRUE;
	  if ((sm&1) == 0) return FALSE;
	  i = 3;
	  while (i*i <= sm )
		{ if ((sm%i) == 0) return FALSE;
		  i++; i++;
		}
	  return TRUE;
 } else if (IS_EVEN(m) || IS_NEG(m)) {
          return FALSE;
 } else {
	  basic_t u[ZSIZEP], minus1[ZSIZEP];
          basic_tc sm = (m[0] == 1 ? m[1] : RADIX);
				/* Upper bound on random number */
          zmont_t mont;
          int i, s;

	  zsadd(m,(long)-1,u);
          s = zremove2(u);	/* Divide  m-1  by power of 2 */
          zmstart(m, &mont);
	  zsub(m, mont.one, minus1);	/* -1 in Montgomery representation */
	  for (i=t; i>0; i--) {
	          basic_t a[ZSIZEP];
		  zintoz(2 + zrandom(sm - 3), a);	   /* 2 <= a < sm */
				        /* (in Montgomery representation) */
/* TBD -- When t > 1, use Lucas test half the time */

		  zmexp(a, u, &mont, a);
		  if (    zcompare(a, mont.one) != 0
		       && zcompare(a, minus1  ) != 0) {	/* Not +1 or -1 yet */
			  int j;
		          boolean_t hit_minus1 = FALSE;
			  for (j = s-1; j > 0 && !hit_minus1; j--)  {
				  zmont(a, a, &mont, a);
				  if (zcompare(a, minus1) == 0) {
				      hit_minus1 = TRUE;
				  }
			  }
			  if (!hit_minus1) return FALSE;
	          }
          } /* for i */
	  return TRUE;	/* Passed t tests -- assumed to be prime */
  } /* else */
} /* zprime */

void zpstart()
{
  if (pshift == -2) {
      long j, jstep, jstart;
      for (j = PRIMBND-1; j >= 0; j--) xxprimes[j] = TRUE;
      for (jstart = 3, jstep = 3;
	   jstart < PRIMBND;
	   jstart += 2*jstep + 2, jstep += 2) {
					/* jstart = (jstep^2 - 3)/2 */

	   if (xxprimes[(jstep-3) >>1]) {	/* If jstep is prime */
		for (j=jstart; j<PRIMBND; j+=jstep) xxprimes[j] = FALSE;
	   }				/* Flag 2*j + 3 as composite */
      }
  }
  pshift = 0;
  pindex = -1;
}

static void zpshift()
{ long j;
  if (pshift == 0) {
	  for (j=PRIMBND-1; j >= 0; j--) interval[j] = xxprimes[j];

  } else {
	  long jstep, jstepsq;
          longc pshift3 = pshift + 3;
	  for (j = PRIMBND-1; j >= 0; j--) interval[j] = TRUE;
	  for (jstepsq = 9, jstep = 3;
	       jstepsq < pshift3 + 2*PRIMBND;
	       jstepsq += 4*jstep + 4, jstep += 2) {
					/* jstepsq = jstep^2 */
	     if (xxprimes[(jstep-3) >> 1]) {
	         long jstart;

		 jstart = 2*jstep - 1 - (pshift3 + jstep - 1) % (2*jstep);
		 jstart = MAX(jstart, jstepsq - pshift3)/2;
		 for (j = jstart; j < PRIMBND; j += jstep) interval[j] = FALSE;
			/* Flag pshift3 + 2*j as composite */
	     } /* if  */
	  } /* for jstep */
  } /*else */
}

long zpnext()	/* returns next prime number, starting at 2 */
		/* restart with pshift = -1 (get 2 first) */
		/* or with zpstart() (get 3 first) */
		/* Primes may grow as large as 4*PRIMBND^2 + 2*PRIMBND */

{ if (pshift < 0) { zpstart(); return(2); }
  if (pindex < 0) { pindex = 0; zpshift(); return(3); }
  while ((++pindex) < PRIMBND ) {
	 if (interval[pindex])  return (pshift+2*pindex+3);
  } /* while */
  pshift += 2*PRIMBND;
  if (pshift > 2*PRIMBND*(2*PRIMBND+1)) { zpstart(); return(2); }
  zpshift();
  pindex = -1;
  while ((++pindex) < PRIMBND )		/* assuming interval non-empty */
	{ if (interval[pindex]) return(pshift+2*pindex+3); }
/*NOTREACHED*/	/* to shut up lint */
}

long isqrt(n)
longc n;
{ long a;
  if (n <= 0) return (0);
  if (n <= 3) return (1);
  if (n <= 8) return (2);
  a = n / 2;
  while (TRUE) {
	longc ndiva = n/a;
	long  newa = (ndiva+a)/2;
	if (newa-ndiva <= 1)
		if (newa*newa <= n) return(newa);
		else return(ndiva);
	a = newa;
	}
/*NOTREACHED*/	/* to shut up lint */
}


void zsqrt(n,r)
basic_tc n[]; basic_t r[];
{ long a[ZZSIZEP], ndiva[ZZSIZEP], diff[ZZSIZEP], i;
  if (IS_ZERO(n)) {
        SET_ZERO(r);
  } else if (IS_NEG(n)) {
	zhalt(ZHALT_ZSQRT_NEGATIVE);
  } else if (n[ZLNGWRD] == 1) {
	zintoz(isqrt(n[1]), r);
  } else	{
	a[0]=(n[0]+1)/2;
	for (i=1;i<a[0];i++) a[i] = 0;
	a[a[0]] = isqrt(n[n[0]]) + 1;
	if ((n[0]&1) ==0) a[a[0]] <<= NBITSH;
	if (a[a[0]]&RADIX) {
		a[a[0]] = 0; a[0] ++; a[a[0]] = 1;
	}
	while (1) {
		zdiv(n,a,ndiva,r);
		zadd(a,ndiva,r);
		i=zsdiv(r,(long)2,r);
		if (zcompare(r,ndiva) <= 0) break;
		zsub(r,ndiva,diff);
		if ((diff[0] == 1) && (diff[1] <= 1)) {
			zmul(r,r,diff);
			if (zcompare(diff,n) > 0) zcopy(ndiva,r);
			break;
			}
		zcopy(r,a);
		}
	}
} /* zsqrt */

