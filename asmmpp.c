/************************************************************************/
/*																		*/
/*		Multi-Precise Arithmetic Package								*/
/*																		*/
/*	Author: Robert D. Silverman											*/
/*		MITRE Corp.														*/
/*		Burlington Rd.													*/
/*		Bedford, Mass.													*/
/*																		*/
/*	Notes:																*/
/*	All multi-precise numbers are stored in radix "RADIX"				*/
/*	They are stored with their length in word 0. The length				*/
/*	includes word 0.													*/
/*	They are stored with most significant digits to the right and		*/
/*	least to the left.													*/
/*																		*/
/*	SIZE(x)  is the size of the multi-precise number x					*/
/*	LASTDIGIT(x)	 is the least significant word						*/
/*	FIRST(x) is the most significant									*/
/*  This version calls (presumed) assembler routines for double length	*/
/*  arithmetic.															*/
/*																		*/
/************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpdim.h"

#define	DEBUG	0

void add(int *a, int *b, int *c),
	 add_single(int *a, int *b, int *c),
	 mult(int *a, int *b, int *c),
	 square(unsigned int *a, unsigned int *c),
	 mpcopy(int *a, int *b),
	 subt(int *a, int *b, int *c),
	 mul_single(int *a, int b, int *c),
     divv(int *a, int *b, int *q, int *r),
	 mpctoi(char *s, int *a),
	 write_num(int *a, char *s),
	 get_num(int *a, char *s),
	 mpsqrt(int *a, int *b),
	 gcd(int *a, int *b, int *c),
	 modinverse(int *a, int *b, int *c),
	 dump_num(int *a, char *s),
	 oldsquare(int *a, int *b),
	 lshift(int *a, int k, int *b),
	 rshift(int *a, int k, int *b),
	 recurs(int a, int b,  int n, int *c),
	 fileit(int *a, FILE *p),
	 mpower(int *a, int *e, int *n, int *r),
	 cpower(int a, int e, int *r),
	 binary(int n, int *bits, int *m),
	 internal(int *a),
	 single_modinv_dp(int a, int *n, int *r);

int spower(int *a, int e, int n),
	mod_mult(int a, int b, int c),
	binary_inv(int a, int n),
	fast_ptest(int *a), 
	ptest(int *a), 
	prp2(int *a), 
	pollardrho(int *a, int *b, int n);

extern void _inline 
					old_div4by2_asm(int *a, int *b, int *r),
					mul_single_prec_asm(int a, int b, int *c),
				    muladd_asm(int a, int b, int c, int *d), 
					divrem_asm(int a, int b, int c, int *d), 
					new_div4by2_asm(int *a, int *b, int *r), 
					square_dble_prec(int *a, int *b),
					square_dble_prec_asm(int *a, int *b);

extern int _inline mod_mult_asm(int a, int b, int c);

/************************************************************************/
/*																		*/
/*						GLOBALS											*/
/*																		*/
/************************************************************************/
 
int ini;
static int MASK[] = { 0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047,
		      4095, 8191, 16383, 32767, 65535, 131071, 262143, 524287,
			  1048575, 2097151, 4194303, 8388607, 16777215, 33554431,
			  67108863, 134217727, 268435455, 536870911, 1073741823 };
 
static int POW2[] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
		      4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288,
			  1048576, 2097152, 4194304, 8388608, 16777216, 33554432,
			  67108864, 134217728, 268435456, 536870912, 1073741824 };
 
static int ONE[MPDIM]  = { 2, 1, (MPDIM-2)*0 };
 
/************************************************************************/
/*																		*/
/*	routine to copy A into B											*/
/*																		*/
/************************************************************************/

void mpcopy(a,b)
int a[],b[];
{   /* start of mpcopy */

register int i;

for (i=0; i<SIZE(a); i++) b[i] = a[i];

}   /* end of mpcopy */



/************************************************************************/
/*																		*/
/*	routine to add A and B	: OUTPUT CAN BE INPUT						*/
/*																		*/
/************************************************************************/
 
void add(a,b,c)
int a[],b[],c[];

{   /* start of add */

register int i,carry,sa;

carry = 0;
sa = SIZE(a);
if (sa == SIZE(b))
	{
	for (i=1; i<sa; i++)
		{
		if ((c[i] = a[i] + b[i] + carry) < RADIX)
		   {
		   carry = 0;
		   }
		else
		   {
		   c[i] -= RADIX;
		   carry = 1;
		   }
		}
	if (carry == 1)
		{
		SIZE(c) = sa + 1;
		FIRST(c) = 1;
		}
	else SIZE(c) = sa;
	}
 

else if (sa < SIZE(b)) add_single(a,b,c);
else add_single(b,a,c);
 
}   /* end of add */

/************************************************************************/
/*																		*/
/*	auxiliary routine for ADD only: assumes SIZE(a) < SIZE(b)			*/
/*	 OUTPUT CAN BE INPUT												*/
/************************************************************************/
 
void add_single(a,b,c)
int a[],b[],c[];
 
{   /* start of add_single */
register int i,carry,as,bs;
 
carry = 0;
as = SIZE(a);
bs = SIZE(b);
for (i=1; i<as; i++)
	{
	if ((c[i] = a[i] + b[i] + carry) < RADIX)
		{
		carry = 0;
		}
	else
		{
		c[i] -= RADIX;
		carry = 1;
		}
	}
 
for (; i< bs; i++)
	{
	if ((c[i] = b[i] + carry) < RADIX) carry = 0;
	else
		{
		c[i] -= RADIX;
		carry = 1;
		}
	}
 
if (carry == 1)
	{
	SIZE(c) = bs + 1;
	FIRST(c) = 1;
	}
else SIZE(c) = bs;

}   /* end of add_single */
 
/************************************************************************/
/*																		*/
/*	routine to compute C = A - B: assumes A >= B						*/
/*	OUTPUT CAN BE INPUT													*/
/************************************************************************/
 
void subt(a,b,c)
int a[],b[],c[];

{   /* start of subt */

register int i,borrow,as,bs;
 
borrow = 0;
as = SIZE(a);
bs = SIZE(b);
for (i=1; i<bs; i++)
	{
	if ((c[i] = a[i] - b[i] - borrow) < 0)
		{
		c[i] += RADIX;
		borrow = 1;
		}
	else borrow = 0;
	}

for (; i<= as; i++)
	{
	if ((c[i] = a[i] - borrow) < 0)
		{
		c[i] += RADIX;
		borrow = 1;
		}
	else borrow = 0;
	}
 
SIZE(c) = as; 
i = as-1;
while (i > 1 && c[i] == 0)
	{
	SIZE(c)--;
	i--;
	}

}   /* end of subt */
 
/************************************************************************/
/*																		*/
/*	routine to do single precision multiply								*/
/*																		*/
/************************************************************************/
 
void mul_single(a,b,c)
int a[],b,c[];
 
{   /* start of mul_single */
register int carry,i,d;
int  x[2];

d = SIZE(a);
carry = 0;
for (i=1; i<d; i++)
	{
	muladd_asm(a[i],b,carry,x);		/* x = a[i]*b + carry		*/
	carry = x[1];
	c[i] = x[0];
	}

if (carry == 0) SIZE(c) = d;
else
	{
	SIZE(c) = d+1;
	FIRST(c) = carry;
	}
 
}   /* end of mul_single */


/************************************************************************/
/*																		*/
/*	routine to do single precision divide								*/
/*																		*/
/************************************************************************/
 
div_single(a,b,result)
unsigned int a[],b,result[];
 
{   /* start of div_single */
register unsigned int i,r,as;
unsigned int  x[2];
 
r = 0;
as = SIZE(a);
for (i=as-1; i>0; i--)
	{
	divrem_asm(r,a[i],b,x);
	result[i] = x[0];
	r = x[1];
	}
SIZE(result) = as;
i = SIZE(result);
while (i > 1 && result[i-1] == 0)
	{
	SIZE(result)--;
	i = SIZE(result);
	}
 
return(r);

}   /* end of div_single */

/************************************************************************/
/*																		*/
/*	routine to multiply C = A*B											*/
/*	OUTPUT CAN'T BE INPUT												*/
/************************************************************************/
 
void mult(a,b,c)
int a[],b[],c[];
 
{   /* start of mult */
register int carry,i,j,k;
int x[2],da,db;
 
da = SIZE(a);
db = SIZE(b);
for (i=1; i< da; i++) c[i] = 0;		/* initialize output		*/
 
for (j=1; j<db; j++)
	{
	carry = 0;
	for (i=1; i<da; i++)
		{
		k = i + j - 1;
		muladd_asm(a[i],b[j],c[k],x);
		x[0] += carry;
		if (x[0] >= RADIX)
		   {
		   x[0] -= RADIX;
		   x[1]++;
		   }
		carry = x[1];
		c[k] = x[0];
		}
	k = da + j - 1;
	c[k] = carry;
	}
 
SIZE(c) = da + db - 1;
if (FIRST(c) == 0) SIZE(c)--;
 
}   /* end of mult */


/************************************************************************/
/*																		*/
/*		routine to perform division										*/
/*																		*/
/************************************************************************/
 
void divv(a,b,q,r)
int a[],b[],q[],r[];
 
{   /* start of divv */
register int i,j,k; 
int qhat,d,m,db,t[2],len;
int x[2],borrow,tmp;
int junk1[MPDIM],junk2[MPDIM],junk3[MPDIM];

/*		if divisor is single precision call div_single			*/

if (SIZE(b) < DOUBLE)
	{
	SIZE(r) = SINGLE;
	r[1] = div_single(a,b[1],q);
	return;
	}

/*		now we do full blown int division (expensive!)			*/ 
 
m = SIZE(a) - SIZE(b) + 2;
db = SIZE(b);
if (m > 1)
	{   /* outer loop */
	d = RADIX/(1 + FIRST(b));
	mul_single(a,d,junk1);
	mul_single(b,d,junk2);		/* normalize input parameters	*/
	junk2[SIZE(junk2)] = 0;		/* pad with high order zero		*/
	len = SIZE(junk1) - 1;
	if (len == SIZE(a)-1) { len++; junk1[len] = 0; }

	for (j=0; j<m-1; j++)
		{						/* compute qhat first			*/
		if (junk1[len] == junk2[db-1])
		   {
		   qhat = RADIX-1;
		   i = junk1[len-1] + junk2[db-1];
		   }
		else
		   {
		   divrem_asm(junk1[len],junk1[len-1],junk2[db-1],t);
		   qhat = t[0];
		   i = t[1];
		   }
		muladd_asm(qhat,junk2[db-2],0,x);

		while (x[1] > i || (x[1] == i && x[0] > junk1[len-2]))
			{
			qhat--;
			i += junk2[db-1];
			muladd_asm(qhat,junk2[db-2],0,x);
			}
		mul_single(junk2,qhat,junk3);	/* multiply and subtract*/

		if (SIZE(junk2) == SIZE(junk3))
			{
			k = SIZE(junk2);
			SIZE(junk3) = k+1;
			junk3[k] = 0;
			}
		borrow = 0;
		k = len - SIZE(junk2) + 1;

		for (i=1; i<SIZE(junk3); i++)	/* subtract				*/
			{
			tmp = junk1[k] - junk3[i] - borrow;
			if (tmp < 0) { borrow = 1; junk1[k] = tmp+RADIX; }
			else         { borrow = 0; junk1[k] = tmp; }
			k++;
			}

		k = m-j-1;						/* test remainder		*/
		if (borrow == 0)  q[k] = qhat; 
		else
			{							/* add back the borrow	*/
			q[k] = qhat-1;
			borrow = 0;
			k = len - SIZE(junk2) + 1;
			for (i=1; i<SIZE(junk2); i++)
			    {
			    tmp = junk1[k] + junk2[i] + borrow;
			    if (tmp < RADIX) { borrow=0; junk1[k]=tmp; }
			    else { borrow=1; junk1[k] = tmp-RADIX; }
			    k++;
			    }
			}
		len--;
		}   /* end loop on j */

	while (db > 2 && junk1[db-1] == 0) db--;
	SIZE(junk1) = db;
	i = div_single(junk1,d,r);			/* unnormalize			*/
	if (q[m-1] == 0) SIZE(q) = m-1;
	else SIZE(q) = m;
	}
else
	{									/* q = 0, r = a			*/
	mpcopy(a,r);
	q[0] = 2;
	q[1] = 0;
	}

}   /* end of div */

 
/************************************************************************/
/*																		*/
/*	routine to convert multi-precision integer to string of				*/
/*	digits of length MPITOC <= size										*/
/*																		*/
/************************************************************************/
 
mpitoc(input,str,size)
int input[];
char str[];
int size;
 
{   /* start of mpitoc */

register int i,j,temp;
unsigned int junk1[MPDIM],d;
static char digits[10] = { '0','1','2','3','4','5','6','7','8','9' };
char tstr[MPCDIM];

mpcopy(input,junk1);

for (i=0; i<size; i++)
	{
	tstr[i] = EOS;
	str[i] = EOS;
	}
i = 0;
while(i < size)
	{
	d = div_single(junk1,10,junk1);
	if (d > 9 || d < 0) { printf("d error %d\n",d); exit(0); }
	tstr[i] = digits[d];
	if ((SIZE(junk1) <= 2 && junk1[1] == 0)) break;
	i++;
	}
 
temp = strlen(tstr);
for (j=0; j<temp; j++) str[j] = tstr[temp-j-1];	/* reverse string	*/
 
return(temp);

}   /* end of mpitoc */

/************************************************************************/
/*																		*/
/*	routine to convert string of digits: input to MP	 				*/
/*	intger output.														*/
/*																		*/
/************************************************************************/
 
void mpctoi(input,output)
int output[];
char input[];

{   /* start of mpctoi */

int d,
    index,
    x[2];

index = 0;
SIZE(x) = 2;
while (input[index] == BLANK) index = index + 1;

SIZE(output) = 2;
for (output[1] = 0; input[index] != EOS; index++)
	{
	d = input[index] - '0';
	mul_single(output,10,output);
	x[1] = d;
	add(output,x,output);
	}
 
}   /* end of mpctoi */

/************************************************************************/
/*																		*/
/*	routine to output an MP integer to the terminal						*/
/*	and to an output string												*/
/*																		*/
/************************************************************************/
 
void write_num(n,numstr)
int n[];
char numstr[];
 
{   /* start of write_num */
register int i;
int nd;
char junkc[MPCDIM];
 
nd = mpitoc(n,junkc,MPCDIM);

for (i=0; i<nd; i++)
	{
	numstr[i] = junkc[i];
	printf("%c",junkc[i]);
	}
numstr[nd] = '\0';
printf("\n");
 
}   /* end of write_num */

/************************************************************************/
/*																		*/
/*	routine to read an MP integer from the terminal						*/
/*																		*/
/************************************************************************/
 
void get_num(number,numstr)
int number[];
char numstr[];
 
{   /* start of get_num */
 
scanf("%s",numstr);					/* read in string of digits	*/
mpctoi(numstr,number);			/* convert to internal format	*/ 
 
}   /* end of get_num */
 
 

/************************************************************************/
/*																		*/
/*	routine to compare two MP integers									*/
/*  returns 1 if a > b,  -1 if a < b, 0 otherwise						*/
/************************************************************************/
 
mpcmp(a,b)
int a[],b[];
 
{    /* start of mpcmp */
register int i;

if (SIZE(a) > SIZE(b)) return(1);
if (SIZE(b) > SIZE(a)) return(-1);
 
for (i=SIZE(a)-1; i > 0; i--)
	{
	if (a[i] > b[i]) return(1);
	if (a[i] < b[i]) return(-1);
	}
 
return(0);

}   /* end of mpcmp */

/************************************************************************/
/*																		*/
/*		routine to return [sqrt(n)]										*/
/*																		*/
/************************************************************************/
 
void mpsqrt(a,b) 
int a[],b[];
 
{   /* start of mpsqrt */

int junkc[MPDIM],junk3[MPDIM];
register int i,da,db;
 
da = SIZE(a);
db = 1 + (da >> 1);
SIZE(b) = db;

for (i=1; i<db; i++) b[i] = 0;

if ((da & 1) == 0) FIRST(b) = 1 + (int)sqrt((double)a[da - 1]);
else FIRST(b) = 1 + (int)(sqrt((double)a[da - 1]*(double)RADIX+(double)a[da - 2]));

if (FIRST(b) >= RADIX)
	{										/* this is rare		*/
	b[db] = 1;
	FIRST(b) -= RADIX;
	SIZE(b) = db + 1;
	}

while(1)									/* loop to infinity	*/
	{
	divv(a,b,junkc,junk3);					/* junkc = a/b		*/
	add(b,junkc,junkc);						/* junkc += b		*/
	da = div_single(junkc,2,junkc);			/* divide by 2		*/
	if (mpcmp(junkc,b) >= 0) break;			/* converged		*/
	mpcopy(junkc,b);
	}
 
}   /* end of mpsqrt */
 
/************************************************************************/
/*																		*/
/*	this routine removes the largest possible power of 2 from N			*/
/*																		*/
/************************************************************************/
 
mpodiv(a)
int a[];
 
{   /* start of mpodiv */
register int i,j,index;
int d,accum;

if (SIZE(a) < 2 && a[1] == 0) return(-1);
i=1;
while (a[i] == 0) i++;
i -= 1;

if (i > 0)
	{
	d = SIZE(a) - i;
	for (j=1; j<=d; j++)
		{
		index = j + i;
		a[j] = a[index];
		}
	SIZE(a) = d;
	}

accum = 1; 
j = 0;
d = a[1];
while ((d & 1) == 0) { j++; d = (d >> 1); accum = (accum << 1); }
if (j > 0) div_single(a,accum,a);

return(MPBITS*i + j);
 
}   /* end of mpodiv */
 
/************************************************************************/
/*																		*/
/*		routine to compute GCD(a,b)										*/
/*																		*/
/************************************************************************/
 
void gcd(a,b,c)
int a[],b[],c[];
 
{   /* start of gcd */
register int i;
int junk1[MPDIM],temp[MPDIM];
 
if (mpcmp(a,b) == 0)					/* if the two are equal	*/
	{
	mpcopy(a,c);
	return;
	}
 
mpcopy(a,c);
mpcopy(b,junk1);
(void) mpodiv(c);
if (SIZE(b) == 2 && b[1] == 0)			/* b is zero			*/
	{
	printf("Attempt to take GCD with 0\n");
	return;
	}

if (SIZE(a) == 2 && a[1] == 0)			/* a is zero			*/
	{
	printf("Attempt to take GCD with 0\n");
	return;
	}

i = mpcmp(c,junk1);
while (i != 0)
	{
	if (i < 0)
		{
		subt(junk1,c,temp);
		mpcopy(temp,junk1);
		(void) mpodiv(junk1);
		}
	else
		{
		subt(c,junk1,temp);
		mpcopy(temp,c);
		(void) mpodiv(c);
		}
	i = mpcmp(c,junk1);
	}
 
}   /* end of gcd */

/************************************************************************/
/*																		*/
/*	routine to compute b^e MOD m for MP integers b,e, and m				*/
/*																		*/
/************************************************************************/
 
void mpower(b,e,m,a)
int b[],e[],m[],a[];
 
{   /* start of mpower */

int junkc[MPDIM],junk4[MPDIM];
int bits[MPBITS+1],ebit;
register int i,n,j,es;

mpcopy(b,a);
es = SIZE(e)-1;
for (n=es; n>0; n--)
	{
	ebit = e[n];
	i = 0;
	for (j=1; j<=MPBITS; j++)
		{
		i++;
		bits[i] = ebit & 1;
		ebit = ebit >> 1;
		if (n == SIZE(e)-1 && ebit == 0) { i--; break; }
		}
	for (; i>0; i--)
		{
		square(a,junkc);
		divv(junkc,m,junk4,a);
		if (bits[i] != 0)
			{
			mult(a,b,junkc);
			divv(junkc,m,junk4,a);
			}
		}
	}
 
}   /* end of mpower */

/************************************************************************/
/*																		*/ 
/*	routine to compute A**N for single precision A and N				*/
/*	and a MP result														*/
/*																		*/
/************************************************************************/
 
void cpower(base,power,result)
int base,
    power,
    result[];
 
{   /* start of cpower */
int temp[MPDIM],
    temp1[MPDIM],
    bits[LOGTWO],
    max;
register int i;
 
INIT(result);
INIT(temp);
binary(power,bits,&max);
SIZE(temp) = 2;
temp[1] = base;
SIZE(result) = 2;
result[1] = base;

for (i=max-1; i>=0; i--)
	{
	if (bits[i] == 1)
		{
		square(result,temp1);
		mul_single(temp1,base,result);
		}
	else
		{
		square(result,temp1);
		mpcopy(temp1,result);
		}
	}

 }  /* end of cpower */


/************************************************************************/
/*																		*/
/*	routine to return binary representation of an integer				*/
/*																		*/
/************************************************************************/
 
void binary(n,bits,max)
int n,bits[],*max;
 
{   /* start of binary */
register int i,p,q,r;
 
for (i=0; i<32; i++) bits[i] = 0;
 
p = n;
while (p > 0)
	{
	q = 1;
	r = 0;
	while (q <= (p >> 1))
		{
		q = q << 1;
		r++;
		}
	bits[r] = 1;
	p = p-q;
	}
 
i = 31;
while (bits[i] != 1) i--;
*max = i;
 
}   /* end of binary */
 

/************************************************************************/
/*																		*/
/*	routine to display a number in internal format on the sceen			*/
/*																		*/
/************************************************************************/
 
void internal(num)
int num[];
 
{   /* start of internal */
int i;
 
for (i=SIZE(num)-1; i>0; i--) printf("%d, ",num[i]);
printf("\n");
 
}   /* end of internal */

/************************************************************************/
/*																		*/
/*	routine to dump current number into a file							*/
/*																		*/
/************************************************************************/
 
void fileit(number,ptr)
int number[];
FILE *ptr;
{
char junk[MPCDIM];
 
(void) mpitoc(number,junk,MPCDIM);
fprintf(ptr,"%s\n",junk);
} 

/************************************************************************/
/*																		*/
/*			prime testing routine; one of several						*/
/*																		*/
/************************************************************************/
 
prime_test(number)
int number[];

{   /* start of prime_test */
 
int temp[MPDIM],temp1[MPDIM],base[MPDIM],q[MPDIM],nminus1[MPDIM],
    quo[MPDIM],k,x,j,i;

x = 101;								/* base for exponentiation	*/
mpcopy(number,temp);					/* make temporary copy		*/
 
SIZE(base) = 2;
FIRST(base) = x;						/* use x as base for test	*/

if ((LASTDIGIT(number) & 1) == 0) return(1);	/* make sure it's odd	*/

temp[1]--;								/* use N-1 as the exponent	*/
mpcopy(temp,nminus1);					/* make copy of N-1			*/

mpower(base,temp,number,temp1);			/* base^(N-1) MOD N			*/

/*	check the residue: if it isn't 1 the number is definitely		*/
/*	composite. Otherwise proceed with the test						*/

if (SIZE(temp1) > 2 || LASTDIGIT(temp1) != 1) return(1);

/*	start by removing largest possible power of 2 from N-1			*/
/*	compute k and q for N = 1 + q*2^k								*/

k = 0;
while (LASTDIGIT(temp) & 1 == 0)
	{
	div_single(temp,2,temp1);
	mpcopy(temp1,temp);
	k++;
	}
mpcopy(temp,q);
 
for (i=0; i<1; i++)
   {
   x = (x*x + 1) % 15511;				/* choose x	'at random'		*/
   INIT(temp1);
   FIRST(base) = x;
   j = 0;
   mpower(base,q,number,temp1);			/* compute x^q MOD N		*/
   while (j < k)						/* strong prp test			*/
	{
	if ( (j == 0 && SIZE(temp1) == 2 && LASTDIGIT(temp1) == 1) ||
	    mpcmp(temp1,nminus1) == 0) goto outer_loop;
	j++;
	square(temp1,temp);
	divv(temp,number,quo,temp1);
	}
   return(1);							/* N is not prime			*/
 
outer_loop:;
   }
 
return(0);								/* 	N is  prime				*/
 
}   /* end of prime_test */
 

/************************************************************************/
/*																		*/
/*	routine to solve x^2 = N mod p, N and p are given					*/
/*	Note: first pass at this routine: use naiive approach				*/
/*																		*/
/************************************************************************/

modsqrt(N,p)
int N,p;
 
{   /* start of modsqrt */
int base[MPDIM];
 
int i;
 
if (N < 0) N += p;
 
if ((p % 4) == 3 && (p > 500))			/* special case 	*/
	{
	SIZE(base) = 2;
	LASTDIGIT(base) = N;
	return(spower(base,(p+1)/4,p));
	}
 
for (i=0; i<p; i++)
	{
	if (i < SQRT_RAD && (i*i % p) == N) return(i);
	if (i > SQRT_RAD && mod_mult_asm(i,i,p) == N) return(i);
	}
 
return(0);		/* no solution		*/
 
}   /* end of modsqrt */
 
/************************************************************************/
/*																		*/
/*	routine to compute base^exp MOD mod for MP integer base				*/
/*	and single precision exp and mod.	 								*/
/*																		*/
/************************************************************************/
 
spower(base,exp,mod)
int base[],exp,mod;
 
{   /* start of spower */

int junk[MPDIM];
int bits[MPBITS+1];
register int i,j,base0,base1;

base1 = div_single(base,mod,junk);
base0 = base1;
 
i = 0;
for (j=1; j<=MPBITS; j++)				/* find binary rep of exp	*/
	{
	i++;
	bits[i] = exp & 1;					/* exp MOD 2				*/
	exp = exp >> 1;
	if (exp == 0) { i--; break; }
	}
for (; i>0; i--)
	{
	if (base1 < SQRT_RAD)
		base1 = base1*base1 % mod;
	else
		base1 = mod_mult_asm(base1,base1,mod);
	if (bits[i] != 0)
		{
		if (base1 < SQRT_RAD && base0 < SQRT_RAD)
			base1 = base0*base1 % mod;
		else base1 = mod_mult_asm(base1,base0,mod);
		}
	}
 
return(base1);
 
}   /* end of spower */
 

 
/************************************************************************/
/*																		*/
/*	routine to return the next prime after N							*/
/*																		*/
/************************************************************************/
 
next_prime(N)
int N;
 
{   /* start of next_prime */
register int p,q,r;
 
if (N == 2) return(3);
if (N == 3) return(5);
 
p = N;
while (1)
	{
	p += 2;
	if ((p % 3) == 0) continue;
	if ((p % 5) == 0) continue;
	q = 5;
	r = (int)sqrt((double)p) + 1;
	while (q <= r)
		{
		q += 2;
		if (p % q == 0) goto nextp;
		q += 4;
		if (p % q == 0) goto nextp;
		}
	return(p);
nextp:;
	}
 
}   /* end of next_prime */

/************************************************************************/
/*																		*/
/*		routine to compute a recursive sequence; naiive					*/
/*																		*/
/************************************************************************/
 
void recurs(s0,s1,term,result)
int result[],
    s0,s1,term;
 
{   /* start of recurs */
register int i;
int temp[MPDIM],temp1[MPDIM];
 
INIT(temp);
INIT(temp1);
SIZE(temp) = 2;
SIZE(temp1) = 2;
temp[1] = s0;
temp1[1] = s1;
 
for (i=2; i<term; i++)
	{
	add(temp,temp1,result);
	mpcopy(temp1,temp);
	mpcopy(result,temp1);
	}
 
}   /* end of recurs */
 

/************************************************************************/
/*																		*/
/*	routine to get single precision remainder							*/
/*	from multi-precise dividend and single divisor						*/
/*																		*/
/************************************************************************/
 
smod(a,b)
int a[],b;
 
{   /* start of single_mod */
register int i,r,as;
int  x[2];
 
r = 0;
as = SIZE(a);
for (i=as-1; i>0; i--)
	{
	divrem_asm(r,a[i],b,x);
	r = x[1];
	}
 
return(r);

}   /* end of single_mod */
 
/************************************************************************/
/*																		*/
/*	routine to shift MP number left by k bits							*/
/*																		*/
/************************************************************************/
 
void lshift(number,k,result)
int number[],result[],k;
 
{   /* start of lshift */
int topk[MPDIM];					/* hold top k bits of each word	*/
int bottom[MPDIM];					/* holds remaining bottom bits	*/
register int i,len,s;
 
len = 30-k;
s = SIZE(number);
for (i=1; i<SIZE(number); i++)
   {
   bottom[i] = number[i] & MASK[len];
   topk[i] = number[i]  >> len;
   }
 
/*		Now see if the size of the number increases					*/
 
result[1] = bottom[1] << k;			/* first element: special case	*/
if (FIRST(number) >= POW2[len])
   {								/* Size increases by 1			*/
   SIZE(result) = s+1;
   for (i=2; i<SIZE(result)-1; i++)
	{
	result[i] = topk[i-1] + (bottom[i] << k);
	}
   result[s] = topk[s-1];			/* last element: special case	*/
   }
else
   {
   SIZE(result) = s;
   for (i=2; i<SIZE(result); i++)
	{
	result[i] = topk[i-1] + (bottom[i] << k);
	}
   }
 }   /* end of lshift */
 
/************************************************************************/
/*																		*/
/*	routine to shift MP number right by k bits							*/
/*																		*/
/************************************************************************/
 
void rshift(number,k,result)
int number[],result[],k;
 
{   /* start of rshift */
 
int topk[MPDIM];				/* hold top k bits of each word	*/
int bottom[MPDIM];				/* and remaining bottom bits	*/
register int i,len,s;
 
 
/*	First split each word into the top k and bottom 30-k bits	*/
 
s = SIZE(number);
len = 30-k;
for (i=1; i<s; i++)
   {
   bottom[i] = number[i] & MASK[k];
   topk[i] = number[i]  >> k;
   }
 
/*		Now see if the size of the number decreases				*/
 
if (topk[s-1] == 0)					/* number reduces in size	*/
   {
   SIZE(result) = s-1;
   for (i=1; i<SIZE(result); i++)
	{
	result[i] = topk[i] + (bottom[i+1] << len);
	}
   }
else
   {
   SIZE(result) = s;
   for (i=1; i<s-1; i++)
	{
	result[i] = topk[i] + (bottom[i+1] << len);
	}
   result[s-1] = topk[s-1];
   }
 
}   /* end of rshift */
 
/************************************************************************/
/*																		*/
/*		Routine to Square an MP number									*/
/*																		*/
/************************************************************************/

void oldsquare(a,c)
int a[],c[];
 
{   /* start of square */
register int carry,i,j,k;
int x[2],da;
 
da = SIZE(a);
for (i=1; i< da; i++) c[i] = 0;		/* initialize output		*/
 
for (j=1; j<da; j++)
	{
	carry = 0;
	for (i=1; i<da; i++)
		{
		k = i + j - 1;
		muladd_asm(a[i],a[j],c[k],x);
		x[0] += carry;
		if (x[0] >= RADIX)
		   {
		   x[0] -= RADIX;
		   x[1]++;
		   }
		carry = x[1];
		c[k] = x[0];
		}
	k = da + j - 1;
	c[k] = carry;
	}
 
SIZE(c) = (da << 1)  - 1;
if (FIRST(c) == 0) SIZE(c)--;
 
}   /* end of square */
 
/************************************************************************/
/*																		*/
/*	routine to compute base^exp MOD mod for  integer base				*/
/*	and single precision exp and mod.	 								*/
/*																		*/
/************************************************************************/
 
fastpow(base,exp,mod)
int base,exp,mod;
 
{   /* start of fastpow */

int bits[MPBITS+1];
register int i,j,base0;

base0 = base;
 
i = 0;
for (j=1; j<=MPBITS; j++)			/* find binary rep of exp	*/
	{
	i++;
	bits[i] = (exp & 1);			/* exp MOD 2				*/
	exp = exp >> 1;
	if (exp == 0) { i--; break; }
	}
for (; i>0; i--)
	{
	if (base < SQRT_RAD)
		base = base*base % mod;
	else
		base = mod_mult_asm(base,base,mod);
	if (bits[i] != 0)
		{
		if (base < SQRT_RAD && base0 < SQRT_RAD)
			base = base0*base % mod;
		else base = mod_mult_asm(base,base0,mod);
		}
	}
 
return(base);
 
}   /* end of fastpow */
 
 
/************************************************************************/
/*																		*/
/*		routine to compute modular inverses for non-prime				*/
/*		moduli via euclidean algorithm									*/
/*																		*/
/************************************************************************/

void modinverse(a,modulus,result)
int a[],modulus[],result[];

{   /* start of modinverse */
 
int ps[MPDIM], ps1[MPDIM], ps2[MPDIM],
    divisor[MPDIM], dividend[MPDIM], rem[MPDIM], 
    q[MPDIM], temp[MPDIM];
register int parity,j;

COPY1(ONE,ps2,j);					/* Use euclidean algorithm	*/
divv(modulus,a,q,rem);				/* initialize the sequence	*/
COPY1(a,dividend,j);
COPY1(rem,divisor,j);
COPY1(q,ps1,j);
parity = 0;
 
while ((SIZE(rem) != SINGLE) || (SIZE(rem) == SINGLE && LASTDIGIT(rem) != 0))
	{								/* while rem != 0			*/
	divv(dividend,divisor,q,rem);	/* q = [dividend/divisor]	*/
	mult(q,ps1,temp);
	add(ps2,temp,ps);				/* ps = q*ps1 + ps2			*/

	COPY1(ps1,ps2,j);				/* Update the sequence		*/
	COPY1(ps,ps1,j);
	COPY1(divisor,dividend,j);
	COPY1(rem,divisor,j);
 
	parity = 1-parity;
	}
 
if (parity == 0) COPY1(ps2,result,j);
else		subt(modulus,ps2,result);
 

}   /* end of modinverse */
 
/************************************************************************/
/*																		*/
/*	routine to output an MP integer to the terminal						*/
/*	and to an output string												*/
/*																		*/
/************************************************************************/
 
void dump_num(n,numstr)
int n[];
char numstr[];
 
{   /* start of write_num */
register long i;
long nd;
char junkc[MPCDIM];
 
nd = mpitoc(n,junkc,MPCDIM);

for (i=0; i<nd; i++)
	{
	numstr[i] = junkc[i];
	printf("%c",junkc[i]);
	}
numstr[nd] = '\0';
 
}   /* end of write_num */
 
/************************************************************************/
/*                                                                      */
/*      Routine to factor a number via pollard-rho                      */
/*                                                                      */
/************************************************************************/

pollardrho(n,fac1,numtrials)
int n[];
int fac1[], numtrials;
 
{   /* start of pollardrho */
int i,j,comp;
int xi[MPDIM],x2i[MPDIM],dif[MPDIM],res[MPDIM];
int temp[MPDIM],temp1[MPDIM];

i = 0;
j = 100;
SIZE(dif) = SINGLE;
FIRST(dif) = 1;
SIZE(xi) = SINGLE;
FIRST(xi) = 1;
SIZE(x2i) = SINGLE;
FIRST(x2i) = 1;

while (1)
   {
   for (; i<j; i++)
        {
        square(xi,temp);
        divv(temp,n,temp1,xi);
        LASTDIGIT(xi)++;
        square(x2i,temp);
        divv(temp,n,temp1,x2i);
        LASTDIGIT(x2i)++;
        square(x2i,temp);
        divv(temp,n,temp1,x2i);
        LASTDIGIT(x2i)++;
        if ((comp = mpcmp(x2i,xi)) == 1) subt(x2i,xi,res);
        else if (comp == -1) subt(xi,x2i,res);
        else return(0);
        mult(res,dif,temp);
        divv(temp,n,temp1,dif);
        if (SIZE(dif) < SINGLE || SIZE(dif) == SINGLE && LASTDIGIT(dif) == 0)
           {
           gcd(res,n,temp);
	   if (SIZE(temp) >= 2 && FIRST(temp) > 1)
		{
		mpcopy(temp,fac1);
		return(1);
		}
           }
        }
   gcd(dif,n,temp);
   if (SIZE(temp) >= SINGLE && FIRST(temp) > 1)
        {
	mpcopy(temp,fac1);
	return(1);
        }
   j = i+120;
   if (j > numtrials) return(0);
   }


}   /* end of pollardrho */


/************************************************************************/
/*                                                                      */
/*      Probable prime routine specifically for double precision        */
/*	actually, this assumes 'number' is at most single-and-a-half		*/
/*	simply tests if number is base 2 Euler psp. Faster than a more		*/
/*	general routine would be.	Good only to 53 bits					*/
/*                                                                      */
/************************************************************************/

prp2(number)
int number[];

{   /* start of prp2 */
int temp[MPDIM], temp1[MPDIM], count, power[MPDIM], bits[128];
int i,ebit;

if (SIZE(number) == SINGLE)					/* EXTREMELY RARE		*/
   {										/* should never happen	*/
   i = fastpow(2,LASTDIGIT(number)-1,LASTDIGIT(number));
   if (i == 1) return(1);
   else return(0);
   }

power[0] = number[0];
power[1] = number[1];
power[2] = number[2];
LASTDIGIT(power)--;
SIZE(temp) = SINGLE;
FIRST(temp) = 2;

ebit = LASTDIGIT(power);				/* convert exponent to binary	*/
count = 0;
for (i=0; i<MPBITS; i++)
   {
   bits[count] = ebit & 1;
   ebit = ebit >> 1;
   count++;
   }
ebit = FIRST(power);
while (ebit > 0)
   {
   bits[count] = ebit & 1;
   ebit = ebit >> 1;
   count++;
   }

for (i=count-2; i>=count-4; i--)	/* first 3 mults need no remainder		*/
   {								/* since they are less than RADIX		*/
   FIRST(temp) = FIRST(temp) * FIRST(temp);
   if (bits[i] == 1) FIRST(temp) = FIRST(temp) << 1;
   }

for (; i>=1; i--)				/* Euler psp is (n-1)/2 for exponent		*/
   {							/* so can ignore last bit					*/
   square_dble_prec_asm(temp,temp1);

   if (bits[i] == 1) mul_single(temp1,2,temp1);	/* shift would be faster	*/

   new_div4by2_asm(temp1,number,temp);
   }

/*	div4by2 is not totally accurate (but it is fast). It might not		*/
/*	return the smallest reduced residue. hence, we do one final			*/
/*	division.															*/

new_div4by2_asm(temp,number,temp1);

/*		number is prime if result is 1 or -1 mod n						*/

if (SIZE(temp1) == SINGLE)
	{
	if (LASTDIGIT(temp1) == 1) return(1);
	if (LASTDIGIT(temp1) == LASTDIGIT(number) - 1) return(1);
	}	
if (FIRST(temp1) == FIRST(number) && LASTDIGIT(temp1) == LASTDIGIT(number) - 1) return(1);
return(0);


}   /* end of prp2 */


/************************************************************************/
/*                                                                      */
/*      Pollard rho for double precision numbers ONLY                   */
/*                                                                      */
/************************************************************************/

pollardrho2(number,fac1,fac2)
int number[], *fac1, *fac2;

{   /* start of pollardrho2 */
int i,j,comp;
int xi[MPDIM],x2i[MPDIM],dif[MPDIM],res[MPDIM];
int temp[MPDIM],temp1[MPDIM];

i = 0;
j = 100;
SIZE(dif) = SINGLE;
FIRST(dif) = 1;
SIZE(xi) = SINGLE;
FIRST(xi) = 1;
SIZE(x2i) = SINGLE;
FIRST(x2i) = 1;

while (1)
   {
   for (; i<j; i++)
	{
	square_dble_prec_asm(xi,temp);
	new_div4by2_asm(temp,number,xi);
	LASTDIGIT(xi)++;
	square_dble_prec_asm(x2i,temp);
	new_div4by2_asm(temp,number,x2i);
	LASTDIGIT(x2i)++;
	square_dble_prec_asm(x2i,temp);
	new_div4by2_asm(temp,number,x2i);
	LASTDIGIT(x2i)++;
	if ((comp = mpcmp(x2i,xi)) == 1) subt(x2i,xi,res);
	else if (comp == -1) subt(xi,x2i,res);
	else return(0);									/* rare			*/
	mult(res,dif,temp);
	new_div4by2_asm(temp,number,dif);
	if (SIZE(dif) < SINGLE || SIZE(dif) == SINGLE && LASTDIGIT(dif) == 0)
	   {
	   gcd(res,number,temp);
	   if (SIZE(temp) != SINGLE || FIRST(temp) == 1) return(0);
	   *fac1 = FIRST(temp);
	   div_single(number,FIRST(temp),temp1);
	   if (SIZE(temp1) > SINGLE) return(0);
	   *fac2 = FIRST(temp1);
	   return(1);
	   }
	}
   gcd(dif,number,temp);
   if (SIZE(temp) == SINGLE && FIRST(temp) > 1)
	  {
	  *fac1 = FIRST(temp);
	  div_single(number,FIRST(temp),temp1);
	  if (SIZE(temp1) > SINGLE) return(0);
	  *fac2 = FIRST(temp1);
	  return(1);
	  }
   j = i+120;
   if (j > 40000) return(0);
   }
}   /* end of pollardrho2 */

/************************************************************************/
/*                                                                      */
/*      Routine to square a multi-precise integer; faster than mult()	*/
/*      see also square_dble_prec										*/
/*																		*/
/************************************************************************/

void square(a,result)
unsigned int a[],result[];

{   /* start of square */
int len, i, j, column, columnp1, ans[2];

len = SIZE(a);
for (i=0; i<=2*len; i++) result[i] = 0;

for (i=1; i<len; i++)
   {
   for (j=i; j<len; j++)
	  {								/* a[i] * a[j] and accumulate	*/
	  column = i + j - 1;
	  columnp1 = i + j;

	  mul_single_prec_asm(a[i],a[j],ans);	/* ans[0] = low half		*/

	  result[column] += ans[0];			/* ans[1] = high half		*/
	  if (result[column] >= RADIX)
		 {
		 result[column] -= RADIX;
		 result[columnp1]++;
		 }
	  result[columnp1] += ans[1];
	  if (result[columnp1] >= RADIX)
		 {
		 result[columnp1] -= RADIX;
		 result[column+2]++;
		 }
	  if (i == j) continue;			/* add twice unless i = j		*/

	  result[column] += ans[0];
	  if (result[column] >= RADIX)
		 {
		 result[column] -= RADIX;
		 result[columnp1]++;
		 }
	  result[columnp1] += ans[1];
	  if (result[columnp1] >= RADIX)
		 {
		 result[columnp1] -= RADIX;
		 result[column+2]++;
		 }
	  }   /* end of j loop */
   }      /* end of i loop */

SIZE(result) = 2*len-1;
if (FIRST(result) == 0) SIZE(result)--;

}   /* end of square */
 
 
/************************************************************************/
/*                                                                      */
/*      Probable prime routine specifically for double precision        */
/*	actually, this assumes 'number' is at most single-and-a-half		*/
/*	simply tests if number is base 29 psp. Faster than a more general	*/
/*	routine would be. Slower than new_prp2								*/
/*                                                                      */
/************************************************************************/

old_prp2(number)
int number[];

{   /* start of old_prp2 */
int temp[MPDIM], temp1[MPDIM], count, power[MPDIM], bits[128];
int i,ebit;

if (SIZE(number) == SINGLE)
   {
   i = fastpow(29,LASTDIGIT(number)-1,LASTDIGIT(number));
   if (i == 1) return(1);
   else return(0);
   }
power[0] = number[0];
power[1] = number[1];
power[2] = number[2];
LASTDIGIT(power)--;
SIZE(temp) = SINGLE;
FIRST(temp) = 29;

ebit = LASTDIGIT(power);				/* convert exponent to binary	*/
count = 0;
for (i=0; i<MPBITS; i++)
   {
   bits[count] = ebit & 1;
   ebit = ebit >> 1;
   count++;
   }
ebit = FIRST(power);
while (ebit > 0)
   {
   bits[count] = ebit & 1;
   ebit = ebit >> 1;
   count++;
   }

if (DEBUG) 
   {
   (void) printf("bits in prp: ");
   for (i=count-2; i>=0; i--) (void) printf("%d",bits[i]);
   (void) printf("\n");
   }

for (i=count-2; i>=0; i--)
   {
   if (DEBUG) { (void) printf("PRP; i = %d %d\n",i,bits[i]); fflush(stdout); }
   if (DEBUG) { (void) printf("temp = "); write_num(temp,numstr); fflush(stdout);}

   square_dble_prec_asm(temp,temp1);

   if (DEBUG) { (void) printf("squared: "); write_num(temp1,numstr); fflush(stdout);}

   if (bits[i] == 1) mul_single(temp1,29,temp1);

   new_div4by2_asm(temp1,number,temp);

   if (DEBUG) { (void) printf("reduced : "); write_num(temp,numstr); fflush(stdout);}
   }

if (SIZE(temp) == SINGLE  && LASTDIGIT(temp) == 1) return(1);
else return(0);


}   /* end of old_prp2 */



/************************************************************************/
/*																		*/
/*	binary modular inverse												*/
/*	only for odd primes and n > 0. 										*/
/*																		*/
/************************************************************************/


int binary_inv(int n, int p)
 
{   /* start of binary_inv */
register int n1,n2, preg, m1, m2;

/*     m1*n == n1 mod p,  m2*n == n2 mod p			*/
/*     n2 is odd, n1 is odd after initialization	*/
 
n1 = n;
n2 = p;
preg = p;
m1 = 1;
m2 = 0;

while (!(n1 & 1)) 
	{
    n1 >>= 1;
	m1 += preg & -(m1 & 1);
//    if (m1 & 1) m1 += preg;
    m1 >>= 1;
	}
if (n1 == 1) return m1;
while (n1 != n2) 
	{
    if (n1 >n2) 
		{
        n1 -= n2;
		m1 -= m2;
        if (m1 < 0) m1 += preg;
        do 
			{
            n1 >>= 1;
			m1 += preg & -(m1 & 1);
//            if (m1 & 1) m1 += preg;
            m1 >>= 1;
			} while (!(n1 & 1));
        if (n1 == 1) return m1;
		} 
	else
		{
        n2 -= n1;
		m2 -= m1;
        if (m2 < 0) m2 += preg;
        do 
			{
            n2 >>= 1;
			m2 += preg & -(m2 & 1);
//            if (m2 & 1) m2 += preg;
            m2 >>= 1;
			} while (!(n2 & 1));
        if (n2 == 1) return m2;
		}
	}
return 0;
}   /* end of binary_inv */



/************************************************************************/
/*																		*/
/*	returns Multiplicative inverse of a mod modulus; a single prec		*/
/*	modulus is double-precision											*/
/*																		*/
/************************************************************************/
 
void single_modinv_dp(a,modulus,result)
int a,modulus[],result[];

{   /* start of single_modinv_dp */
 
int ps[MPDIM], ps1[MPDIM], ps2[MPDIM],
    divisor[MPDIM], dividend[MPDIM], rem[MPDIM], 
    q[MPDIM], temp[MPDIM], a1[MPDIM];
register int parity,j;

SIZE(a1) = SINGLE;
LASTDIGIT(a1) = a;

COPY1(ONE,ps2,j);					/* Use euclidean algorithm	*/
divv(modulus,a1,q,rem);				/* initialize the sequence	*/
COPY1(a1,dividend,j);
COPY1(rem,divisor,j);
COPY1(q,ps1,j);
parity = 0;
 
while ((SIZE(rem) != SINGLE) || (SIZE(rem) == SINGLE && LASTDIGIT(rem) != 0))
	{								/* while rem != 0			*/
	divv(dividend,divisor,q,rem);	/* q = [dividend/divisor]	*/
	mult(q,ps1,temp);
	add(ps2,temp,ps);				/* ps = q*ps1 + ps2			*/

	COPY1(ps1,ps2,j);				/* Update the sequence		*/
	COPY1(ps,ps1,j);
	COPY1(divisor,dividend,j);
	COPY1(rem,divisor,j);
 
	parity = 1-parity;
	}
 
if (parity == 0) COPY1(ps2,result,j);
else		subt(modulus,ps2,result);

}   /* end of single_modinv_dp */


/************************************************************************/
/*																		*/
/*	Fast modular square root routine: assumes sqrt exists				*/
/*																		*/
/*	Essentially uses Shanks method										*/
/*																		*/
/************************************************************************/
 
msqrt(u,p)
int u,p;
 
{   /* start of msqrt */
 
register int i,v,w,r,e,q,h,k,y;
int yp,z;
static int two[] = {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,
					32768,65536,131072,262144,524288,1048576};
 
 
if (u < 0) u += p;			/* make sure u is positive	*/
 
if (p < 1000)				/* special case for small p and	*/
   {						/* small prime powers			*/
   for (i=1; i<p; i++) if (i*i % p == u) return(i);
   return(0);				/* didn't find one		*/
   }
 
e = 1;
q = (p-1)/2;
for (i=0; i<20; i++)		/* find q,e such that u = q*2^e + 1	*/
   {
   if ((q & 1) == 0) {q = q >> 1; e++;}
   }
 
 
v = fastpow(u,(q+1)/2,p);		/* u^(q+1)/2 and u^q+1  mod p	*/
w = fastpow(u,q,p); 
 
if (w == 1) return(v);
 
for (i=2; i<=50; i++)			/* find a non-residue of p	*/
	{
	z = fastpow(i,q,p);
	if (fastpow(z,two[e-1],p) == p-1) break;
	}
 
y = z;
r = e;
while(1)				/* find the order of u		*/
   {
   h = w;
   for (k=1; k<=e; k++) if ((h = mod_mult_asm(h,h,p)) == 1) break;
 
   yp = y;
 
   for (i=1; i<r-k; i++) yp = mod_mult_asm(yp,yp,p);
 
   v = mod_mult_asm(v,yp,p);
   yp = mod_mult_asm(yp,yp,p);
   y = yp;
   r = k;
   w = mod_mult_asm(w,yp,p);
 
   if (w == 1) return(v);
   }
 
}   /* end of msqrt */


/************************************************************************/
/*																		*/
/*			prime testing routine; one of several; robust but slow		*/
/*																		*/
/************************************************************************/
 
ptest(number)
int number[];

{   /* start of ptest */
 
int temp[MPDIM],temp1[MPDIM],base[MPDIM];
										/* base for exponentiation	*/
mpcopy(number,temp);					/* make temporary copy		*/
 
SIZE(base) = 2;
FIRST(base) = 101;						/* use 101 as base for test	*/

if ((LASTDIGIT(number) & 1) == 0) return(1);	/* make sure it's odd	*/

temp[1]--;								/* use N-1 as the exponent	*/

mpower(base,temp,number,temp1);			/* base^(N-1) MOD N			*/

/*	check the residue: if it isn't 1 the number is definitely		*/
/*	composite. Otherwise proceed with the test						*/

if (SIZE(temp1) > 2 || LASTDIGIT(temp1) != 1) return(0);
else return(1);							/* prime					*/

}   /* end of ptest */

/************************************************************************/
/*																		*/
/*	faster ptest; if N < 2^53 call prp2;								*/
/*																		*/
/************************************************************************/

fast_ptest(number)
int number[];
 
{   /* start of fast_ptest */
#define EBASE  23
int junkc[MPDIM],junk4[MPDIM],e[MPDIM];
int bits[MPBITS+1],ebit, a[MPDIM];
register int i,n,j,es;

if (FIRST(number) < (1<<22)) return(prp2(number));

SIZE(e) = SIZE(number);
LASTDIGIT(e) = LASTDIGIT(number);
FIRST(e) = FIRST(number);
LASTDIGIT(e)--;

SIZE(a) = SINGLE;
LASTDIGIT(a) = EBASE;

es = SIZE(e)-1;
for (n=es; n>0; n--)
	{
	ebit = e[n];
	i = 0;
	for (j=1; j<=MPBITS; j++)
		{
		i++;
		bits[i] = ebit & 1;
		ebit = ebit >> 1;
		if (n == SIZE(e)-1 && ebit == 0) { i--; break; }
		}
	for (; i>0; i--)
		{
		square_dble_prec_asm(a,junkc);
		divv(junkc,number,junk4,a);
		if (bits[i] != 0)
			{
			mul_single(a,EBASE,junkc);
			divv(junkc,number,junk4,a);
			}
		}
	}

if (LASTDIGIT(a) == 1) return(1);
else return(0);

}   /* end of fast_ptest */
