#include "dimods.h"

/*
	This file contains additions to the basic.c subroutines.
	The file is #included when basic.c is compiled,
	so both groups appear on basic.o.
*/

#if 0
void dtosz(d, a) doublec d; CONSTP(SMPNUM) a;
		/* Convert double precision to nearest SMPNUM. */
{
  if (d >= 0) {
    a->sign = 1;
    dtoz(d, a->abs);
  } else {
    a->sign = -1;
    dtoz(-d, a->abs); 
  }
  if (IS_ZERO(a->abs)) {
    SSET_ZERO(a);
  }
}
#endif

void dtoz(d, a) doublec d; basic_t a[];
	/* Convert nonnegative double precision to nearest MPNUM. */
{
  static doublec dradix = (double)RADIX;
  static doublec dradix_inv = 1.0/(double)RADIX;
  register double da = fabs(d) + 0.5;
  if (da < 1.0) {
    SET_ZERO(a);
  } else {
    register int lng = 1;
    while (da >= dradix) {
      lng++;
      da *= dradix_inv;
    }
    a[ZLNGWRD] = (d > 0  ? lng : -lng);;
    while (lng > 0) {
      register basic_tc digit = (basic_t)da;

      a[lng] = digit;
      da = dradix * (da - (double)digit);
      lng--;
    }
  } /* else */
}

int maxpw2(n) basic_tc n;
/*	Return largest power of 2 dividing a nonzero integer. */
{
    register ubasic_t nn = (ubasic_t)n;
    register int pw2 = 0;
    static char lowpw2[8] = {3, 0, 1, 0, 2, 0, 1, 0};

    if (nn == 0) zhalt(ZHALT_MAXPW2_ZERO);

    while ((nn & 15) == 0) {pw2 += 4; nn >>= 4;}
    pw2 += lowpw2[nn & 7];
    return pw2;
}


basic_t seinv(n, p) basic_tc n, p;
{ 
  if (n > 0 && n < p && (p & 1) != 0) {
    return sinv(n, p);	/* See basic.c */
  } else { 
    register long q=n, r, nq, nr, u=1, w, nw, par=1;
    w = p/n;
    r = p - w*n; 
    while (r != 0) {	/* u*n == par*q   and   w*n == -par*r  mod p */
      nq = q/r;         nr = q - nq*r;
      nw = nq*w + u;	u = w;
       w = nw;		q = r;
       r = nr;		par = -par;
    }
    if (q != 1)  zhalt(ZHALT_SEINV_NOT_COPRIME);

    return (par == 1) ? u : p-u;
  }
}

basic_t sqrtp(nn, p)
basic_tc nn, p;
{
/*
		This routine returns the square root of nn modulo p.
		p should be a prime, and nn should be a quadratic
		residue modulo p.
*/
    basic_tc pm1half = (p-1)/2;
    basic_t n;		/* nn mod p */
    basic_t k;

    n = nn % p;
    if (n < 0) n += p;

    if (p == 2 || n == 0) return n;

    if ((p & 3) == 3) return zsexp (n, (p+1)/4, p);

    for (k = 1; (1L << (k/2)) < p; k++) {
				/* somewhat arbitrary upper bound */
	basic_t xcoef = 1;
	basic_t ccoef = k;
	basic_t ibit;
/*
		    Compute (x+k)^( (p-1)/2 ) mod (x^2 - n)

		    where k is semi-random and x is an indeterminate.

		    We hope GCD ( (x+k)**( (p-1)/2 ) - 1, x**2 - n)
		    has degree 1, in which case its root is our answer.

*/
	ibit = 1L;
	while (2*ibit <= pm1half) ibit *= 2;

 	while ((ibit >>= 1) != 0) {
	    basic_t itemp;
/*
		Square (xcoef*x + ccoef)
*/
	    itemp =   zsmulmod(xcoef, xcoef, p);
	    xcoef = 2*zsmulmod(xcoef, ccoef, p);
	    if (xcoef >= p) xcoef -= p;

	    ccoef = zsmulmod(ccoef, ccoef, p) + zsmulmod(itemp, n, p);
	    if (ccoef >= p) ccoef -= p;
/*
		Possibly multiply by (x + k).
*/
	    if ((pm1half & ibit) != 0) {
		itemp = xcoef;
		xcoef = ccoef + zsmulmod(xcoef, k, p);
		if (xcoef >= p) xcoef -= p;

		ccoef = zsmulmod(ccoef, k, p) + zsmulmod(itemp, n, p);
		if (ccoef >= p) ccoef -= p;
	    }

/*
			We have (x+k)** FLOOR(pm1half/ibit) ==
				 xcoef*x + ccoef mod (x**2 - n) mod p

			We hope the right side is -1 for one
			value of x and +1 for the other value.
			This can happen only if ccoef = 0 since p is odd.
			Then the equations

				  xcoef*x + ccoef == 1 mod p
				- xcoef*x + ccoef == -1 mod p
				              x*x == n mod p

			imply ccoef = 0 (since p is odd) and
			x == n*xcoef mod p.
*/
	    if (ccoef == 0) {
		basic_tc ians = zsmulmod(n, xcoef, p);
		if (zsmulmod(ians, ians, p) == n) return ians;
	    }
	} /* while */
    } /* for k */

    fprintf(stderr, "sqrtp - Can\'t get square root of %ld mod %ld\n", nn, p);
    zhalt(ZHALT_SQRTP_UNSUCCESSFUL);
    return 0;		/* avoid compiler warning */
}

#if 0
void szadd(a, b, c) CONSTP(SMPNUMC) a; CONSTP(SMPNUMC) b;
		    CONSTP(SMPNUM) c;
{
  if (a->sign + b->sign != 0) {
    c->sign = a->sign | b->sign;      /* If both are alike, keep sign */
				      /* If either is zero, use sign of other */
    zadd(a->abs, b->abs, c->abs);
  } else {
    longc zc = zcompare(a->abs, b->abs);
    if (zc > 0) {
      c->sign = a->sign;
      zsub(a->abs, b->abs, c->abs);
    } else if (zc < 0) {
      c->sign = b->sign;
      zsub(b->abs, a->abs, c->abs);
    } else {
      SSET_ZERO(c);
    }
  }
} /* szadd */


void szcopy(a, b) CONSTP(SMPNUMC) a; CONSTP(SMPNUM) b;
{
  b->sign = a->sign;
  zcopy(a->abs, b->abs);
}


void szmul(a, b, c) CONSTP(SMPNUMC) a; CONSTP(SMPNUMC) b;
		    CONSTP(SMPNUM) c;
{
  c->sign = a->sign * b->sign;
  zmul(a->abs, b->abs, c->abs);
} /* szmul */


basic_t szsdiv(a, sngl, quot) CONSTP(SMPNUMC) a; basic_tc sngl; 
			      CONSTP(SMPNUM) quot;
	/* Divide by single precision number (possibly negative). */
	/* Remainder has same sign as denominator. */
{
 
  register basic_tc idiv = labs(sngl);
  register sign_t qsign = a->sign;
  register basic_t irem = zsdiv(a->abs, idiv, quot->abs);
  			/* First divide by ABS(sngl) */
  if (qsign < 0) {
    irem = -irem;
    if (irem < 0) {
      irem += idiv;
      zsadd(quot->abs, 1L, quot->abs); /* Increment ABS(quotient) */
    }
  }
  if (sngl < 0) {
    qsign = -qsign;
    irem = -irem;	/* Negate quotient, remainder if divisor is negative */
  }
  quot->sign = qsign;
  return irem;
}  




void szsmul(a, sngl, c) CONSTP(SMPNUMC) a; basic_tc sngl; CONSTP(SMPNUM) c;
{		/* c = a * sngl.  sngl may be negative */
  if (sngl > 0) {
    c->sign = a->sign;
    zsmul(a->abs, sngl, c->abs);
  } else if (sngl < 0) {
    c->sign = -(a->sign);
    zsmul(a->abs, -sngl, c->abs);
  } else {
    SSET_ZERO(c);
  }
}


char* szsout(a, str) CONSTP(SMPNUMC) a; CONSTP(char) str; 
{
  if (a->sign < 0) {
    str[0] = '-';
    zsout(a->abs, str+1);
  } else {
    zsout(a->abs, str);
  } 
  return str;
}


void szsread(str, snum) CONSTP(charc) str; CONSTP(SMPNUM) snum;
{		/* Read a number from a string.  No error detection. */
  snum->sign = zsread(str, snum->abs);
  snum->abs[ZLNGWRD] = LABS(snum->abs[ZLNGWRD]);
}


void szsub(a, b, c) CONSTP(SMPNUMC) a; CONSTP(SMPNUMC) b;
		    CONSTP(SMPNUM) c;
{
  if  (a->sign == b->sign) {		/* same sign */
    longc zc = zcompare(a->abs, b->abs);
    if (zc > 0) {
      c->sign = a->sign;
      zsub(a->abs, b->abs, c->abs);
    } else if (zc < 0) {
      c->sign = -(a->sign);
      zsub(b->abs, a->abs, c->abs);
    } else {
      SSET_ZERO(c);
    }
  } else {
    c->sign = a->sign | (-b->sign);  
    zadd(a->abs, b->abs, c->abs);
  } 
} /* szsub */

double sztod(a) CONSTP(SMPNUMC) a;
{
  register double d = ztod(a->abs);
  if (a->sign < 0) d = -d;
  return d;
}
#endif

sign_t zditoz(di, mp) DBLINTC di; MPNUM mp;
{				/* Convert DBLINT to MP */

  if (di == (DBLINT)0) {
    SET_ZERO(mp);
    return 0;
  } else {
    register DBLINT dileft = (di < 0 ? -di : di);
    register int lng;

    lng = 0;
    while (dileft >= (DBLINT)RADIX) {
      register basic_tc lowbits = dimods(dileft, RADIX);

      mp[++lng] = lowbits;
      dileft = (dileft - (DBLINT)lowbits)/(DBLINT)RADIX;
    }
    mp[++lng] = dileft;
    if (di > (DBLINT)0) {
      mp[ZLNGWRD] = lng;
      return 1;
    } else {
      mp[ZLNGWRD] = -lng;
      return -1;
    }
  }
} /* zditoz */


void zdivx(a, b, c)	/* c = a/b, error if remainder nonzero */
  basic_tc a[], b[]; basic_t c[];
{
  MPNUM rem;
  zdiv(a, b, c, rem);
  if (!IS_ZERO(rem)) {
    NUMSTR ostr;
    fprintf(stderr, "zdivx - Nonzero remainder\n");
    fprintf(stderr, "Numerator = %s\n", zsout(a, ostr)); 
    fprintf(stderr, "Denominator = %s\n", zsout(b, ostr)); 
    fprintf(stderr, "Quotient = %s\n", zsout(c, ostr)); 
    fprintf(stderr, "Remainder = %s\n", zsout(rem, ostr)); 
    zhalt(ZHALT_ZDIVX_NONZERO_REMAINDER);
  }
}



void zmexp(base, expon, mont, result)
        basic_tc* base;		/* base */
        basic_tc* expon;	/* exponent */
CONSTP( zmont_tc) mont;	        /* modulus information */
        basic_t*result;	        /* output */
/*
		Compute exponentials using zmont.  Actually sets

		answer = base^expon / RADIX^(length*(expon-1))
		      mod modulus.

		Output can overlap inputs.

		The program pre-computes base^2, base^3, base^5, and base^7.
   	        Then it uses what is essentially a radix-8 algorithm.  
                For large values of expon, it takes about
		LOG2(expon) squarings and LOG2(expon)/4 general
		multiplications.
*/

#define MULTIPLY(a, b, c) zmont(a, b, mont, c)
{
    int lex = expon[0];		/* Remaining exponent length */
    int carry = expon[lex];
    if (carry == 0) {
        zcopy(mont->one, result);
    } else {
        MPNUM pw2, pw3, pw5, pw7;  /* base^2, base^3, base^5, base^7 */
	MPNUM answer;		   /* Final answer */
        int  lg;
	boolean_t have7;	/* True if pw7 is defined */
/*
		Get four leading bits of exponent.  However,
		if expon < 8, set carry to the entire exponent.
*/
	lg = 0;
	while ((carry >> lg) > 1) {lg++;}	/* lg = base 2 log of carry */

	lg -= 3;
	if (lg < 0) {
	    lex--;
	    if (lex == 0) {
		lex = 1;
		lg = 0;
	    } else {
		carry = (carry << (-lg)) | (expon[lex] >> (NBITS+lg));
		lg += NBITS;
	    }
	} else {
	    carry >>= lg;
	}
#if 0
	printf("carry, lg = %ld %d\n", carry, lg);
#endif

/*
	pw2, pw3, pw5, pw7 hold base^2, base^3, base^5, base^7 respectively.
*/

        have7 = FALSE;
	if (carry >= 2) { 
	    MULTIPLY(base, base, pw2);		/* base^2 */
	    if (carry >= 3) {
	        MULTIPLY(base, pw2, pw3);       /* base^3 */
	        if (carry >= 5) {
	            MULTIPLY(pw3, pw2, pw5);    /* base^5 */
	            if (carry == 7 || carry == 9 || carry > 11) {
		        have7 = TRUE;		/* Means pw7 defined */
	    		MULTIPLY(pw5, pw2, pw7);
		    }	
		}
	    }	
	}

/*
		Either initialize answer = base^carry and set carry = 0,
		or initialize answer = base^(carry-1) and set carry = 1,
		or initialize answer = base^(carry-2) and set carry = 2.
*/
	switch(carry) {
	default:
	    fprintf(stderr, "zmexp- carry out of range - 1 %d\n", carry);
	    zhalt(ZHALT_ZMEXP_INVALID_CARRY_1);
	case 1:
	    carry = 0;
	    zcopy(base, answer);
	    break;
	case 2:
	    carry = 0;
	    zcopy(pw2, answer);
	    break;
	case 3:
	case 4:
	    carry -= 3;
	    zcopy(pw3, answer);
	    break;
	case 5:
	case 6:
	    carry -= 5;
	    zcopy(pw5, answer);
	    break;
	case 7:
	case 9:
	    carry -= 7;
	    zcopy(pw7, answer);
	    break;
	case 8:
	    carry = 0;
	    MULTIPLY(pw3, pw5, answer);
	    break;
	case 10:
	case 11:
	    carry -= 10;
	    MULTIPLY(pw5, pw5, answer);
	    break;
	case 12:
	case 13:
	    carry -= 12;
	    MULTIPLY(pw5, pw7,  answer);
	    break;
	case 14:
	case 15:
	     carry -= 14;
	     MULTIPLY(pw7, pw7, answer);
	     break;
	}

	while(TRUE) {
/*
		At this point, the desired result is

		  (answer*(base^carry))^(2^(lg + NBITS*lex))
		* base^(bottom bits of exponent)

		where 0 <= carry <= 3
*/

	    lg--;
	    if (lg < 0) {
		lex--;
		lg = NBITS - 1;
		if (lex == 0) break;
	    }

	    carry = (2*carry) | ((expon[lex]>>lg) & 1);

/*
                       carry   New answer               New carry

                         0     answer^2                    0
                         1     answer^2                    1
                         2     answer^2                    2
                         3     answer^2                    3
                         4     (answer*pw2)^2              0
                         5     (answer^2)*pw5              0
                         6     (answer*pw3)^2              0
			 7     (answer^2*pw7               0
*/
	    switch(carry) {
	    default:
		fprintf(stderr,
			"zmexp - carry out of range - 2 %d\n", carry);
	        zhalt(ZHALT_ZMEXP_INVALID_CARRY_2);
	    case 0:
	    case 1:
	    case 2:
	    case 3:			/* answer = answer^2 */
		MULTIPLY(answer, answer, answer);
		break;
	    case 4:			/* answer = (answer*pw2)^2 */
		carry = 0;
		MULTIPLY(answer, pw2, answer); 
		MULTIPLY(answer, answer, answer);
		break;
	    case 5:			/* answer = (answer^2)*pw5 */
		carry = 0;
		MULTIPLY(answer, answer, answer);
		MULTIPLY(answer, pw5, answer);
		break;
	    case 6:			/* answer = (answer*pw3)^2 */
		carry = 0;
		MULTIPLY(answer, pw3, answer);
		MULTIPLY(answer, answer, answer);
		break;
	    case 7:			/* answer = (answer^2)*pw7 */
		carry = 0;
		if (!have7) {		/* set pw7 = pw5*pw2 if needed */
		    have7 = TRUE;
		    MULTIPLY(pw5, pw2, pw7);
		}
		MULTIPLY(answer, answer, answer);
		MULTIPLY(answer, pw7, answer);
		break;
	    }
	}		/* while(TRUE) */

/*		Final processing.  Multiply answer by base^carry.	*/

	switch(carry) {
	default:
	    fprintf(stderr, "zmexp - carry out of range - 3 %d\n", carry);
	    zhalt(ZHALT_ZMEXP_INVALID_CARRY_3);
	case 0:
	    zcopy(answer, result);
	    break;
	case 1:
	    MULTIPLY(answer, base, result);
	    break;
	case 2:
	    MULTIPLY(answer, pw2, result);
	    break;
	case 3:
	    MULTIPLY(answer, pw3, result);
	    break;
	}
    }			/* lex == 0 */
#undef MULTIPLY
}

#ifndef ZOGCD_PRINTS
#define ZOGCD_PRINTS 0
#endif
boolean_t zogcd(n1, n2, gcd) basic_tc n1[], n2[]; basic_t gcd[];
    /* gcd = odd GCD(n1, n2), where n1, n2 not both even .*/
    /* Returns TRUE if GCD > 1, else FALSE. */
    /* Algorithm has special code for operands of lengths 1 and 2. */

{
    MPNUM a, b, quot;
#if ZOGCD_PRINTS
    NUMSTR o1, o2;
#endif

    zcopyabs(n1, a);
    zcopyabs(n2, b);
    if (IS_EVEN(a) && IS_EVEN(b)) zhalt(ZHALT_ZOGCD_BOTHEVEN);	
    if (IS_ZERO(n1) || IS_ZERO(n2)) {
	zadd(a, b, gcd);
	return !IS_ONE(gcd);	/* One operand is zero, other is the GCD */
    }
    if (IS_EVEN(a)) zremove2(a);
    if (IS_EVEN(b)) zremove2(b);

    while (a[0] > 2 || b[0] > 2) {	/* Until both are double precision */
					/* a and b are odd and positive */
	register intc alng = a[0], blng = b[0];
	register basic_tc alead = a[alng], blead = b[blng];
	register basic_tc sum4tst = (a[1] + b[1]) & 2;
			/* Zero if a + b is divisible by 4 */
#if ZOGCD_PRINTS
	printf("a = %s,    b = %s\n", zsout(a, o1), zsout(b, o2));
#endif
	if (alng > blng) {
	    if (alng > blng + 1 || alead >= (blead >> (NBITS-5))) {
		zdiv(a, b, quot, a);
		if (IS_ZERO(a)) {
		    zcopy(b, gcd); return !IS_ONE(gcd);
		}
	    } else if (sum4tst == 0) {
		zadd(a, b, a);
	    } else {
		zsub(a, b, a);
	    }
	    if (IS_EVEN(a)) zremove2(a);
	} else if (blng > alng) {
	    if (blng > alng + 1 || blead >= (alead >> (NBITS-5))) {
		zdiv(b, a, quot, b);
		if (IS_ZERO(b)) {
		    zcopy(a, gcd); return !IS_ONE(gcd);
		}
	    } else if (sum4tst == 0) {
		zadd(b, a, b);
	    } else {	
		zsub(b, a, b);
	    }
	    if (IS_EVEN(b)) zremove2(b);
	} else {			/* Equal lengths */
	    register basic_t zc = alead - blead;
	    if (zc == 0) zc = zcompare(a, b);
	    if (zc < 0) {	/* a < b */
		if (sum4tst == 0) {
		    zadd(b, a, b);
		} else {
		    zsub(b, a, b);
		}
		zremove2(b);
	    } else if (zc > 0) { /* a > b */
		if (sum4tst == 0) {
		    zadd(a, b, a);
		} else {
		    zsub(a, b, a);
		}
		zremove2(a);
	    } else {		/* a == b, so gcd = a */
		zcopy(a, gcd);
		return TRUE;	
	    }
	}
    }	/* while */
    {
	register basic_t a1 = a[1], b1 = b[1];
	if (a[0] > 1 || b[0] > 1) {	/* Double precision case */
	    register basic_t a2 = (a[0] > 1 ? a[2] : 0);
	    register basic_t b2 = (b[0] > 1 ? b[2] : 0);

	    while (a2 + b2 > 0) {
		register basic_t zc = a2 - b2;

#if ZOGCD_PRINTS
		printf("a = %ld %ld, b = %ld %ld\n", a2, a1, b2, b1);
#endif
		if (zc == 0) {
		    zc = a1 - b1;	
		    if (zc == 0) {
			gcd[0] = 2; gcd[1] = a1; gcd[2] = a2;
			return TRUE;
		    }	
		}
		if (zc > 0) {
		    a2 -= b2;
		    a1 -= b1;
		    if (a1 < 0) {a1 += RADIX; a2--;}
		    do {
			a1 = (a1 >> 1) + ((a2 & 1) << (NBITS-1));
			a2 >>= 1;
		    } while ((a1 & 1) == 0);
		} else {
		    b2 -= a2;
		    b1 -= a1;
		    if (b1 < 0) {b1 += RADIX; b2--;}
		    do {
			b1 = (b1 >> 1) + ((b2 & 1) << (NBITS-1));
			b2 >>= 1;
		    } while ((b1 & 1) == 0);
		} /*if */
	    } /* while a2 + b2 */
	}			/* End double precision code */
	while (a1 != b1) {	/* Single precision GCD loop */
#if ZOGCD_PRINTS
	    printf("a = %ld, b = %ld\n", a1, b1);
#endif
	    if (a1 > b1) {
		a1 -= b1;
		do {
		    a1 >>= 1;
		} while ((a1 & 1) == 0);
	    } else {
		b1 -= a1;
		do {
		    b1 >>= 1;
		} while ((b1 & 1) == 0);
	    }
	}  /* while */
	gcd[0] = 1;
	gcd[1] = a1;	
	return a1 > 1;
    }
}


boolean_t zegcd(n1, n2, gcd) basic_tc n1[], n2[]; basic_t gcd[];
	/* gcd = GCD(n1, n2) */
{

#if 0
    NUMSTR o1, o2;
    printf("zegcd inputs: %s %s\n", zsout(n1, o1), zsout(n2, o2));
#endif
     
    if (IS_ZERO(n1)) {
      zcopyabs(n2, gcd);
    } else if (IS_ZERO(n2)) {
      zcopyabs(n1, gcd);
    } else {
      MPNUM m1, m2;
      int pw1, pw2;

      zcopyabs(n1, m1);   pw1 = zremove2(m1);
      zcopyabs(n2, m2);   pw2 = zremove2(m2);

      pw2 = MIN(pw1, pw2);
      zogcd(m1, m2, gcd);
      while (pw2 > 0) {
	zadd(gcd, gcd, gcd);
	pw2--;
      }
    }
    return !IS_ONE(gcd);
}


int zremove2(n) basic_t n[];
	/* Divide n (which must be nonzero) by as large
	   a power of 2 as possible, return exponent of that power. */
{
    longc slngn = n[0];
    register long lngn = LABS(slngn);
    register int i;
    int pwr2 = 0;

    if (IS_ZERO(n)) zhalt(ZHALT_ZREMOVE2_ZERO); 

    while (n[1] == 0) {
	lngn--;
        pwr2 += NBITS;
	for (i = 1; i <= lngn; i++) {
	    n[i] = n[i+1];		/* Divide by RADIX */
	}
    }
    if ((n[1] & 1) == 0) {
	register basic_tc mask = (n[1] - 1) & ~n[1];
	register intc shft = maxpw2(mask+1);
				/* Power of 2 to divide by */
	register basic_t carry = 0;

	pwr2 += shft;
	for (i = lngn; i > 0; i--) {
	    register basic_t newni = (carry << (NBITS-shft)) | (n[i] >> shft);
	    carry = n[i] & mask;
	    n[i] = newni;
	}
	if (n[lngn] == 0) lngn--;	
    }
    n[0] = (slngn >= 0 ? lngn : -lngn);
    if (IS_EVEN(n)) zhalt(ZHALT_ZREMOVE2_INTERNALERR);
    return pwr2;
}


char *zsout_help(a, buf)
basic_t a[];		/* Modified  (and nonnegative initially) */
char buf [];		/* Output */
{
  register int i, nc = 0;
  register basic_t irem;
#if RADIX < 1000000000
  #error - "zsout_help -- RADIX too small"
#endif

/*
	Do divides by 10^9 until quotient is single precision.
*/
  while (a[ZLNGWRD] > 1) {
    irem = zsdiv(a, 1000000000L, a);
    for (i = 0; i < 8; i++) {
      buf[nc++] = (char)(irem % 10);
      irem /= 10;
    } 
    assert (irem >= 0 && irem < 10);
    buf[nc++] = (char)irem;
  } 

/*
	Convert leading digits to decimal.
*/
  irem = (IS_ZERO(a) ? 0 : a[1]);
  while (irem >= 10) {
    buf[nc++] = (char)(irem % 10);
    irem /= 10;
  }
  buf[nc++] = (char)irem;

  for (i = 0; 2*i <= nc-1; i++) {   /* Reverse digits, add '0' to each */
    charc ch = buf[i] + '0';
    buf[i] = buf[nc-1-i] + '0';
    buf[nc-1-i] = ch;
  }
  buf[nc] = '\0';
  return buf;
}


char *zsout(a,buf) basic_tc a[]; CONSTP(char) buf;
{
  MPNUM ca;
  zcopyabs(a, ca);
  if (IS_NEG(a)) {
    buf[0] = '-';
    zsout_help(ca, buf+1);
  } else {
    zsout_help(ca, buf);
  }
  return buf;
}


sign_t zsread(string, coef) CONSTP(charc) string; MPNUM coef;
   /* Read decimal integer.  No check for integer overflow. */
   /* Return -1 if negative, 0 if zero, +1 if positive, ZSREADERR if error. */
{
    sign_t sign = 1;	/* +1 */
    charc *numstr = string;
    if (*numstr == '-') {
	numstr++;
	sign = -1;
    } else if (*numstr == '+') {
	numstr++;
    }

    SET_ZERO(coef);
    while (*numstr != '\0'){
	basic_tc digit = *numstr++ - '0';
	if (digit < 0 || digit > 9) {
	   return ZSREADERR;	/* Invalid digit */
	} else if (coef[ZLNGWRD] == 1 && coef[1] < RADIX/10) {
	   coef[1] = 10*coef[1] + digit;
	} else {
	   zsmul(coef, (basic_t)10, coef);
	   zsadd(coef, digit,       coef);
	}
    }
    if (sign < 0) coef[ZLNGWRD] = -coef[ZLNGWRD];
    return IS_ZERO(coef) ? 0 : sign;
}


void zsubabs(a, b, c) basic_tc a[], b[]; basic_t c[];
/* c = abs(a-b), output can be input */
{
    if (zcompare(a,b) >= 0) {
	zsub(a, b, c);
    } else {
	zsub(b, a, c);
    }
}


double ztod(mp) basic_tc mp[];
/* Convert multiple precision number to floating point.
   No check is made for exceeding largest allowable FP number.
*/
{
    register double d = 0.0;
    register int i;
    longc slng = mp[ZLNGWRD];

    for (i = LABS(slng); i > 0; i--) {
	d = d*(double)RADIX + (double)mp[i];
    }
    return (slng >= 0 ? d : -d);
}

