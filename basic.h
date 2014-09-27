/*
	This file contains declarations for
	programs using the basic.c multiple
	precision package.
*/
#ifndef BASIC_H
#define BASIC_H 


#include "nfstypes.h"


/* NBITS must be even */
#define NBITS 30
#define NBITSH (NBITS>>1)
#define RADIX (1L<<NBITS)
#define RADIXM (RADIX-1)
#define RADIXH (1L<<(NBITS-1))
#define RADIXROOT (1L<<NBITSH)
#define RADIXROOTM (RADIXROOT-1)
#define ZLNGWRD 0
#define ZSIZEP 23
#define ZZSIZEP (2*ZSIZEP-1)
#define	ZDECSIZE (ZSIZEP*NBITS*301/500 + 1)
		/* 301/500 = 2*log(2)/log(10) */

/*
  zhalt  error codes 0-99 are reserved to basic.c.
	 error codes 100-199 are reserved to this include file and basic2.c
	 error codes 200 and above are for application program.
*/
#define ZHALT_ZMSTART_EVEN_MODULUS 3L
#define ZHALT_ZMSTART_TOO_LONG 4L
#define ZHALT_ZRANDOM_UNINITIALIZED 5L
#define ZHALT_SINV_INVALID_ARGUMENTS 7L
#define ZHALT_SODINV_NOT_COPRIME 8L

#define ZHALT_ZADDPOS_BAD_LENGTHS  9L
#define ZHALT_ZSUBPOS_BAD_LENGTHS 10L
#define ZHALT_ZSUBPOS_NEGATIVE 11L
 /* 12L unused */
#define ZHALT_ZSQRT_NEGATIVE 13L
#define ZHALT_ZDIV_ZERO_DIVISOR 14L
#define ZHALT_ZMOD_ZERO_DIVISOR 15L

#define ZHALT_ZREMOVE2_ZERO 101L
#define ZHALT_ZREMOVE2_INTERNALERR 102L
#define ZHALT_ZOGCD_BOTHEVEN 103L
#define ZHALT_MAXPW2_ZERO 104L
#define ZHALT_SEINV_NOT_COPRIME 105L
#define ZHALT_ZDIVX_NONZERO_REMAINDER 106L
#define ZHALT_ZMEXP_INVALID_CARRY_1  107L
#define ZHALT_ZMEXP_INVALID_CARRY_2  108L
#define ZHALT_ZMEXP_INVALID_CARRY_3  109L
#define ZHALT_SQRTP_UNSUCCESSFUL     110L


#define IS_EVEN(x) (x[ZLNGWRD] == 0 || (x[1] & 1) == 0)
#define IS_ONE(x) (x[ZLNGWRD] == 1 && x[1] == 1)
#define IS_ZERO(x) (x[ZLNGWRD] == 0)
#define IS_NEG(x) (x[ZLNGWRD] < 0)
#define SET_ZERO(x) x[ZLNGWRD] = 0
#define SET_ONE(x)  x[1] = x[ZLNGWRD] = 1
#define SSET_ONE(x) do {(x)->sign = 1; (x)->abs[0] = 1; (x)->abs[1] = 1;}\
		      while(FALSE)
#define SSET_ZERO(x) do {(x)->sign = 0; (x)->abs[0] = 1; (x)->abs[1] = 0;}\
		      while(FALSE)

typedef signed long basic_t;	/* Digit in multiple precision number */
typedef const basic_t basic_tc;	/* Constant basic_t */
typedef unsigned long ubasic_t;
typedef const ubasic_t ubasic_tc;

typedef basic_t MPNUM[ZSIZEP];
typedef basic_tc MPNUMC[ZSIZEP];

typedef int sign_t;		/* Sign of multiple precision number */
				/* -1, 0, +1, or ZSREADERR */
typedef const sign_t sign_tc;
#define ZSREADERR 2

typedef struct {		/* Data for Montgomery multiplicaiton (zmont) */
	ubasic_t zminv;		/* -1/n[1] mod RADIX */
        basic_t n[ZSIZEP];	/* Modulus passes to zmstart */
        basic_t one[ZSIZEP]; 	/* Constant +1 in Montgomery form */
        basic_t tom[ZSIZEP];	/* zmont by this to convert to Montgomery form*/
      } zmont_t;

typedef const zmont_t zmont_tc;

typedef char NUMSTR[ZDECSIZE];


#define MPINPUT basic_tc*
#define MPOUTPUT basic_t*
#define MPMODIFIED basic_t*

        /* -------- Start of routines on basic.c file -------- */

void fzread  ARGS((FILE*, MPOUTPUT));                     /* Read from a file */
void zadd    ARGS((MPINPUT, MPINPUT, MPOUTPUT));                  /* Addition */
void zaddmod ARGS((MPINPUT, MPINPUT, MPINPUT, MPOUTPUT)); /* Modular addition */
void zcopy   ARGS((MPINPUT, MPOUTPUT));                /* Copy 1st arg to 2nd */
void zcopyabs ARGS((MPINPUT, MPOUTPUT));          /* Copy abs(1st arg) to 2nd */
void zdiv    ARGS((MPINPUT, MPINPUT, MPOUTPUT, MPOUTPUT)); /* Division, Q & R */
void zdiv_help ARGS((MPINPUT, longc, MPINPUT, longc, MPOUTPUT, MPOUTPUT)); 
						          /* zdiv/zmod helper */
void zexp    ARGS((MPINPUT, MPINPUT, MPINPUT, MPOUTPUT));   /* Modular expon. */
void zhalt   ARGS((basic_tc));                      	     /* Error handler */
void zintoz  ARGS((basic_tc, MPOUTPUT));        /* Convert integer to MP form */
void zmstart ARGS((MPINPUT, CONSTP(zmont_t)));   /* Montgomery initialization */
void zmont   ARGS((MPINPUT, MPINPUT, CONSTP(zmont_tc), MPOUTPUT));
						  /* Montgomery multplication */
void zmod    ARGS((MPINPUT, MPINPUT, MPOUTPUT));              /* Remaindering */
void zmul    ARGS((MPINPUT, MPINPUT, MPOUTPUT));            /* Multiplication */
void zmulmod ARGS((MPINPUT, MPINPUT, MPINPUT, MPOUTPUT));    /* Modular mult. */
void zpstart ARGS((NOARGS));                    /* Initialize prime generator */
void ztom    ARGS((MPINPUT, CONSTP(zmont_tc), MPOUTPUT));   
						     /* To Montgomery repres. */
void zmtoz   ARGS((MPINPUT, CONSTP(zmont_tc), MPOUTPUT)); 
						   /* From Montgomery repres. */
void zread   ARGS((MPOUTPUT));                    /* Read from standard input */
void zrstart ARGS((basic_tc));                /* Random number initialization */
void zsadd   ARGS((MPINPUT, basic_tc, MPOUTPUT));        /* Add integer to MP */
void zsmul   ARGS((MPINPUT, basic_tc, MPOUTPUT));   /* Multiply MP by integer */
void zsqrt   ARGS((MPINPUT, MPOUTPUT));                        /* Square root */
void zstart  ARGS((NOARGS));                                /* Initialization */
void zsub    ARGS((MPINPUT, MPINPUT, MPOUTPUT));               /* Subtraction */
void zsubmod ARGS((MPINPUT, MPINPUT, MPINPUT, MPOUTPUT));     /* Modular sub. */

basic_t fzwrite  ARGS((FILE*, MPINPUT, basic_tc, charc*, charc*));
                                                             /* Write to file */
basic_t isqrt    ARGS((basic_tc));           /* Truncated integer square root */
basic_t sinv     ARGS((basic_tc, basic_tc)); /* Invert 1st arg modulo odd 2nd */
basic_t zcompabs ARGS((MPINPUT, MPINPUT));  /* Compare absolute values
					               -- output -1, 0, or +1 */
basic_t zcompare ARGS((MPINPUT, MPINPUT));  /* Compare -- output -1, 0, or +1 */
basic_t zinv     ARGS((MPINPUT, MPINPUT, MPOUTPUT));      /* Inversion or gcd */
basic_t zpnext   ARGS((NOARGS));                    /* Next prime in sequence */
boolean_t zprime ARGS((MPINPUT,intc));               /* Probable prime tester */
basic_t zrandom  ARGS((basic_tc));                 /* Random number generator */
basic_t zsdiv    ARGS((MPINPUT, basic_tc, MPOUTPUT)); /* Divide MP by integer */
basic_t zsexp    ARGS((basic_tc, basic_tc, basic_tc));      /* Modular expon. */
basic_t zsmulmod ARGS((basic_tc, basic_tc, basic_tc));       /* Modular mult. */
basic_t ztoint   ARGS((MPINPUT));                    /* Convert MP to integer */
basic_t zwrite   ARGS((MPINPUT));                 /* Write to standard output */


        /* -------- End of routines on basic.c file -------- */

        /* -------- Following routines defined in basic2.c. ------- */


void      dtoz     ARGS ((doublec, MPOUTPUT));	     /* Round double to MPNUM */
int       maxpw2   ARGS((basic_tc));       /* Power of 2 dividing nonzero int */
sign_t    zditoz	ARGS((DBLINTC, MPOUTPUT));    /* Convert DBLINT to mp */
basic_t   seinv    ARGS((basic_tc, basic_tc));  /* sinv allowing even 2nd arg */
basic_t   sqrtp    ARGS ((basic_tc, basic_tc));	/* Square root modula a prime */

void      zdivx    ARGS ((MPINPUT, MPINPUT, MPOUTPUT)); 
				  /* Divide, issue error if remainder nonzero */
boolean_t zegcd    ARGS ((MPINPUT, MPINPUT, MPOUTPUT));   /* Even or odd  GCD */

void zmexp ARGS ((MPINPUT, MPINPUT, CONSTP(zmont_tc), MPOUTPUT));
                                                /* Exponentiation using zmont */
void zmlucas ARGS ((MPINPUT, MPINPUT, CONSTP(zmont_tc), MPOUTPUT));
				/* Lucas function using zmont */

boolean_t zogcd    ARGS ((MPINPUT, MPINPUT, MPOUTPUT));            /* Odd GCD */
int       zremove2 ARGS ((MPMODIFIED)); /* Remove powers of 2 dividing number */
char*     zsout_help  ARGS ((MPMODIFIED, char*));         /* Helper for zsout */
char*     zsout    ARGS ((MPINPUT, CONSTP(char)));      /* Write MP to string */
sign_t    zsread   ARGS ((CONSTP(charc), MPOUTPUT));    /* Read from a string */
void      zsubabs  ARGS ((MPINPUT, MPINPUT, MPOUTPUT));  /* Abs value of diff */
double    ztod     ARGS ((MPINPUT));                  /* Convert MP to double */

        /* -------- End of routines on basic2.c file -------- */

	/* Following routines  defined in assembly language routines */

void      PRDPPA   ARGS ((ubasic_tc, ubasic_tc, ubasic_t[2]));
					/* Add product to third argument. */
					/* Arguments must be nonnegative. */

#endif /* BASIC_H */
